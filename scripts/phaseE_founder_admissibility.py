#!/usr/bin/env python3
"""
Phase-E — Founder-Conditioned Developmental Admissibility (LARRY)

Goal:
  Test whether clone-level founder purity preferentially coincides with clones that
  have low temporal entropy (i.e., temporally coherent developmental progression),
  using a fixed, nonparametric score.

Score definition (no fit, no knobs):
  For each clone c:
    H_c = Shannon entropy of Time point within clone (bits)
    w_c = 1 - H_c / log2(n_times)   (clipped to [0,1])
    purity_c = max_f P(founder=f | clone=c)
  Phase-E score:
    S = mean_c ( w_c * purity_c )

Null:
  Shuffle founder labels within time strata (time label = Time point code) and
  recompute S using the same w_c (computed from the chosen time labels).

Permutation p-value:
  p = (1 + #{S_null >= S_obs}) / (1 + perms)

Decision (one-sided):
  GO if (viability_ok) AND (S_obs > S_null_95) AND (p <= alpha)

Optimization notes:
  - Clone structure is fixed across permutations.
  - Time entropy per clone is computed once for the selected time labels (possibly controlled).
  - Clone purity per permutation is computed via a single np.bincount over (clone,founder).

Control modes:
  - none
      Standard Phase-E.

  - shuffle_time_within_clone
      SANITY / invariance check: permute time labels within each clone.
      Preserves each clone's time histogram => H_c unchanged => w_c unchanged.
      This should NOT "collapse" the score; it's not a destructive control for this metric.

  - shuffle_founder_within_clone
      SANITY / invariance check: permute founder labels within each clone.
      Preserves each clone's founder histogram => purity_c unchanged.
      This should NOT "collapse" the score; it's not a destructive control for this metric.

  - shuffle_time_across_clones
      DESTRUCTIVE control: permute time labels globally across cells.
      Breaks clone↔time coupling and typically reduces score magnitude.
      Does **not necessarily eliminate detectability** if founder structure remains.

  - shuffle_founder_across_clones
      TRUE destructive control: permute founder labels globally across cells.
      Breaks clone↔founder coupling and should eliminate detectability if Phase-E holds.

  - permute_w_across_clones
      DIAGNOSTIC control: keep clone time structure intact (so w_c is meaningful),
      but break alignment between w_c and purity_c by permuting w across clones.
      This should reduce the score if the effect relies on alignment (not just marginals).

  - resample_time_by_clone_from_global
      DESTRUCTIVE control: for each clone, resample time labels IID from the *global* time
      distribution (preserves clone sizes + global time marginal; destroys clone↔time coupling).

Inputs (Phase-1 joined TSV):
  Required columns:
    - clone_id
    - clone_valid
    - Starting population          (founder)
    - Time point                  (time label)

Output (out_dir/):
  - phaseE_summary.txt
  - clone_time_entropy.tsv
  - null_distribution.tsv
  - clone_terms_table.tsv
  - founder_majority_summary.tsv
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Tuple

import numpy as np
import pandas as pd


def _entropy_bits_from_counts(counts: np.ndarray) -> float:
    """Counts -> Shannon entropy in bits. counts is 1D nonnegative."""
    total = float(counts.sum())
    if total <= 0.0:
        return 0.0
    p = counts[counts > 0].astype(np.float64) / total
    return float(-(p * np.log2(p)).sum())


def _factorize_str(s: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
    codes, uniques = pd.factorize(s.astype(str), sort=True)
    return codes.astype(np.int32), np.asarray(uniques, dtype=str)


@dataclass
class Encoded:
    clone: np.ndarray
    founder: np.ndarray
    time: np.ndarray
    n_clones: int
    n_founders: int
    n_times: int
    clone_names: np.ndarray
    founder_names: np.ndarray
    time_names: np.ndarray


def load_larry_encoded(path_tsv: str) -> Encoded:
    df = pd.read_csv(path_tsv, sep="\t")
    needed = ["clone_id", "clone_valid", "Starting population", "Time point"]
    for c in needed:
        if c not in df.columns:
            raise RuntimeError(f"Missing required column: {c}")

    # Phase-E: only valid clones
    df = df[df["clone_valid"] == True].copy()  # noqa: E712
    if len(df) == 0:
        raise RuntimeError("No rows after clone_valid filter.")

    clone_codes, clone_uniques = _factorize_str(df["clone_id"])
    founder_codes, founder_uniques = _factorize_str(df["Starting population"])
    time_codes, time_uniques = _factorize_str(df["Time point"])

    n_clones = int(clone_uniques.size)
    n_founders = int(founder_uniques.size)
    n_times = int(time_uniques.size)

    if n_founders < 2:
        raise RuntimeError("Phase-E viability NO-GO: fewer than two founders after filtering.")

    return Encoded(
        clone=clone_codes,
        founder=founder_codes,
        time=time_codes,
        n_clones=n_clones,
        n_founders=n_founders,
        n_times=n_times,
        clone_names=clone_uniques,
        founder_names=founder_uniques,
        time_names=time_uniques,
    )


def compute_time_entropy_per_clone(enc: Encoded, time_labels: np.ndarray) -> np.ndarray:
    """
    Compute H(Time | clone) in bits, per clone, using a single bincount over (clone,time).
    Returns: time_entropy_bits[clone]
    """
    if time_labels.shape != enc.time.shape:
        raise ValueError("time_labels must have same shape as enc.time")

    nC, nT = enc.n_clones, enc.n_times
    key = enc.clone.astype(np.int64) * nT + time_labels.astype(np.int64)
    ct = np.bincount(key, minlength=nC * nT).reshape(nC, nT)

    out = np.zeros(nC, dtype=np.float64)
    for c in range(nC):
        out[c] = _entropy_bits_from_counts(ct[c])
    return out


def compute_founder_purity_per_clone(enc: Encoded, founder_labels: np.ndarray) -> np.ndarray:
    """
    Given founder_labels per cell, compute purity per clone:
      purity[clone] = max_f count(clone,f) / count(clone)
    Using a single bincount over (clone,founder).
    """
    if founder_labels.shape != enc.founder.shape:
        raise ValueError("founder_labels must have same shape as enc.founder")

    nC, nF = enc.n_clones, enc.n_founders
    key = enc.clone.astype(np.int64) * nF + founder_labels.astype(np.int64)
    cf = np.bincount(key, minlength=nC * nF).reshape(nC, nF)

    clone_sizes = cf.sum(axis=1).astype(np.float64)
    denom = np.maximum(clone_sizes, 1.0)
    purity = cf.max(axis=1).astype(np.float64) / denom
    return purity


def shuffle_founders_within_time(
    founder: np.ndarray, time_labels: np.ndarray, n_times: int, seed: int
) -> np.ndarray:
    """
    Null generator: shuffle founder labels within each time stratum.
    Returns a new founder array of same length.
    """
    rng = np.random.default_rng(seed)
    f = founder.copy()
    for t in range(int(n_times)):
        idx = np.flatnonzero(time_labels == t)
        if idx.size <= 1:
            continue
        f[idx] = rng.permutation(f[idx])
    return f


def shuffle_time_within_clone(clone: np.ndarray, time: np.ndarray, seed: int) -> np.ndarray:
    """
    SANITY control: shuffle time labels within each clone.

    Preserves per-clone time histogram -> preserves per-clone entropy -> preserves w.
    """
    rng = np.random.default_rng(seed)
    out = time.copy()
    for c in np.unique(clone):
        idx = np.flatnonzero(clone == c)
        if idx.size <= 1:
            continue
        out[idx] = rng.permutation(out[idx])
    return out


def shuffle_founder_within_clone(clone: np.ndarray, founder: np.ndarray, seed: int) -> np.ndarray:
    """
    SANITY control: shuffle founder labels within each clone.

    Preserves per-clone founder histogram -> preserves per-clone purity.
    """
    rng = np.random.default_rng(seed)
    out = founder.copy()
    for c in np.unique(clone):
        idx = np.flatnonzero(clone == c)
        if idx.size <= 1:
            continue
        out[idx] = rng.permutation(out[idx])
    return out


def shuffle_time_across_clones(time: np.ndarray, seed: int) -> np.ndarray:
    """DESTRUCTIVE control: permute time labels globally across cells."""
    rng = np.random.default_rng(seed)
    return rng.permutation(time)


def shuffle_founder_across_clones(founder: np.ndarray, seed: int) -> np.ndarray:
    """DESTRUCTIVE control: permute founder labels globally across cells."""
    rng = np.random.default_rng(seed)
    return rng.permutation(founder)


def resample_time_by_clone_from_global(enc: Encoded, seed: int) -> np.ndarray:
    """
    DESTRUCTIVE control:
      For each clone, resample time labels IID from the global time distribution.
    Preserves:
      - clone sizes
      - global time marginal (in expectation; exactly per-cell sampling)
    Destroys:
      - clone↔time coupling
    """
    rng = np.random.default_rng(seed)
    global_time_counts = np.bincount(enc.time, minlength=enc.n_times).astype(np.float64)
    global_time_probs = global_time_counts / global_time_counts.sum()

    out = enc.time.copy()
    for c in np.unique(enc.clone):
        idx = np.flatnonzero(enc.clone == c)
        if idx.size == 0:
            continue
        out[idx] = rng.choice(enc.n_times, size=idx.size, p=global_time_probs)
    return out


def _compute_weights_from_entropy(enc: Encoded, time_entropy_bits: np.ndarray) -> np.ndarray:
    """
    No-knob normalization:
      w = 1 - H / log2(n_times), clipped to [0,1]
    """
    H_max = float(np.log2(max(enc.n_times, 2)))  # >= 1.0
    w = 1.0 - (time_entropy_bits / H_max)
    w = np.clip(w, 0.0, 1.0)
    return w.astype(np.float64)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_tsv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument(
        "--control_mode",
        choices=[
            "none",
            "shuffle_time_within_clone",
            "shuffle_founder_within_clone",
            "shuffle_time_across_clones",
            "shuffle_founder_across_clones",
            "resample_time_by_clone_from_global",
            "permute_w_across_clones",
        ],
        default="none",
        help=(
            "Controls for Phase-E. "
            "SANITY: shuffle_time_within_clone, shuffle_founder_within_clone. "
            "DESTRUCTIVE: shuffle_time_across_clones, shuffle_founder_across_clones, "
            "resample_time_by_clone_from_global. "
            "DIAGNOSTIC: permute_w_across_clones."
        ),
    )
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    enc = load_larry_encoded(args.input_tsv)

    # Viability gate (minimal): enough clones to estimate something meaningful
    viability_ok = enc.n_clones >= 100

    base_seed = int(args.seed)

    # Apply control (once) to define "observed" labels used in scoring and null construction.
    if args.control_mode == "none":
        time_obs = enc.time
        founder_obs = enc.founder
        control_desc = "none"
    elif args.control_mode == "shuffle_time_within_clone":
        time_obs = shuffle_time_within_clone(enc.clone, enc.time, seed=base_seed + 777)
        founder_obs = enc.founder
        control_desc = "SANITY: shuffle time within clone (time histograms preserved)"
    elif args.control_mode == "shuffle_founder_within_clone":
        time_obs = enc.time
        founder_obs = shuffle_founder_within_clone(enc.clone, enc.founder, seed=base_seed + 777)
        control_desc = "SANITY: shuffle founder within clone (founder histograms preserved)"
    elif args.control_mode == "shuffle_time_across_clones":
        time_obs = shuffle_time_across_clones(enc.time, seed=base_seed + 777)
        founder_obs = enc.founder
        control_desc = "DESTRUCTIVE: shuffle time across clones (clone↔time broken)"
    elif args.control_mode == "shuffle_founder_across_clones":
        time_obs = enc.time
        founder_obs = shuffle_founder_across_clones(enc.founder, seed=base_seed + 777)
        control_desc = "DESTRUCTIVE: shuffle founder across clones (clone↔founder broken)"
    elif args.control_mode == "resample_time_by_clone_from_global":
        time_obs = resample_time_by_clone_from_global(enc, seed=base_seed + 777)
        founder_obs = enc.founder
        control_desc = "DESTRUCTIVE: resample time per clone from global time distribution"
    elif args.control_mode == "permute_w_across_clones":
        time_obs = enc.time
        founder_obs = enc.founder
        control_desc = "DIAGNOSTIC: permute w across clones (w↔purity alignment broken)"
    else:
        raise RuntimeError(f"Unknown control_mode: {args.control_mode}")

    # Fixed per-clone time entropy and weights (given time_obs)
    time_entropy = compute_time_entropy_per_clone(enc, time_obs)
    w = _compute_weights_from_entropy(enc, time_entropy)

    # If alignment-control is selected, permute w across clones (keep marginals, break alignment)
    if args.control_mode == "permute_w_across_clones":
        rng = np.random.default_rng(base_seed + 888)
        w = rng.permutation(w)

    # Observed founder purity and score (use founder_obs, not enc.founder)
    purity_obs = compute_founder_purity_per_clone(enc, founder_obs)
    S_obs = float(np.mean(w * purity_obs))

    # --- Interpretability add-on (no knobs, no gating) ---
    nC, nF = enc.n_clones, enc.n_founders
    key_obs = enc.clone.astype(np.int64) * nF + founder_obs.astype(np.int64)
    cf_obs = np.bincount(key_obs, minlength=nC * nF).reshape(nC, nF)

    maj_founder_code = cf_obs.argmax(axis=1).astype(np.int32)
    maj_founder_name = enc.founder_names[maj_founder_code]

    clone_sizes = cf_obs.sum(axis=1).astype(np.float64)
    maj_share = cf_obs.max(axis=1).astype(np.float64) / np.maximum(clone_sizes, 1.0)

    df_clone = pd.DataFrame(
        {
            "clone_id": enc.clone_names,
            "majority_founder": maj_founder_name,
            "clone_size": clone_sizes.astype(np.int32),
            "time_entropy_bits": time_entropy,
            "weight_w": w,
            "founder_purity_obs": purity_obs,
            "majority_share": maj_share,
            "term_w_times_purity": w * purity_obs,
        }
    )

    summ = (
        df_clone.groupby("majority_founder", sort=True)
        .agg(
            n_clones=("clone_id", "count"),
            median_clone_size=("clone_size", "median"),
            mean_clone_size=("clone_size", "mean"),
            median_time_entropy_bits=("time_entropy_bits", "median"),
            mean_time_entropy_bits=("time_entropy_bits", "mean"),
            median_founder_purity=("founder_purity_obs", "median"),
            mean_founder_purity=("founder_purity_obs", "mean"),
            median_term=("term_w_times_purity", "median"),
            mean_term=("term_w_times_purity", "mean"),
        )
        .reset_index()
    )

    df_clone.to_csv(os.path.join(args.out_dir, "clone_terms_table.tsv"), sep="\t", index=False)
    summ.to_csv(os.path.join(args.out_dir, "founder_majority_summary.tsv"), sep="\t", index=False)

    # Null: shuffle founders within time strata (respect time_obs)
    B = int(args.perms)
    S_null = np.zeros(B, dtype=np.float64)

    for b in range(B):
        f_sh = shuffle_founders_within_time(
            founder_obs, time_obs, enc.n_times, seed=base_seed + 1 + b
        )
        purity_b = compute_founder_purity_per_clone(enc, f_sh)
        S_null[b] = float(np.mean(w * purity_b))

        if (b + 1) % max(1, B // 10) == 0:
            print(f"  perm {b+1:>4}/{B}  S_null={S_null[b]:.6f}")

    S_med = float(np.median(S_null))
    S_95 = float(np.quantile(S_null, 0.95))
    p_S = float((1 + np.sum(S_null >= S_obs)) / (1 + B))  # +1 correction

    alpha = float(args.alpha)
    decision = "GO" if (viability_ok and (S_obs > S_95) and (p_S <= alpha)) else "NO-GO (detectability)"

    print("=== Phase-E summary ===")
    print(f"n_clones={enc.n_clones}")
    print(f"n_founders={enc.n_founders}")
    print(f"viability_ok={'true' if viability_ok else 'false'}")
    print(f"control_mode={args.control_mode} ({control_desc})")
    print(f"S_obs={S_obs}")
    print(f"S_null_median={S_med}")
    print(f"S_null_95={S_95}")
    print(f"p_S={p_S}")
    print(f"Decision={decision}")

    # Outputs
    pd.DataFrame(
        {
            "clone_id": enc.clone_names,
            "time_entropy_bits": time_entropy,
            "weight_w": w,
            "founder_purity_obs": purity_obs,
        }
    ).to_csv(os.path.join(args.out_dir, "clone_time_entropy.tsv"), sep="\t", index=False)

    pd.DataFrame({"S_null": S_null}).to_csv(
        os.path.join(args.out_dir, "null_distribution.tsv"), sep="\t", index=False
    )

    with open(os.path.join(args.out_dir, "phaseE_summary.txt"), "w", encoding="utf-8") as f:
        f.write("=== Phase-E summary ===\n")
        f.write(f"n_clones={enc.n_clones}\n")
        f.write(f"n_founders={enc.n_founders}\n")
        f.write(f"viability_ok={'true' if viability_ok else 'false'}\n")
        f.write(f"control_mode={args.control_mode} ({control_desc})\n")
        f.write(f"S_obs={S_obs}\n")
        f.write(f"S_null_median={S_med}\n")
        f.write(f"S_null_95={S_95}\n")
        f.write(f"p_S={p_S}\n")
        f.write(f"Decision={decision}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

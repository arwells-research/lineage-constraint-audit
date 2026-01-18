#!/usr/bin/env python3
"""
phaseD_founder_dispersion.py — Phase-D: Founder dispersion / transition admissibility (LARRY)

Question:
  Do different founders occupy distinct admissible developmental "transition regions,"
  beyond merely differing in total contribution?

LARRY binding (from phase1_joined.tsv header):
  - unit: clone_id (clone_valid == true)
  - founder: Starting population
  - developmental axis: Time point (categorical)
  - dispersion metric (primary): Shannon entropy (bits) of Time point within clone

Test statistic (overall):
  - For each founder f with n_clones >= min_units_per_founder:
        m_f = median(entropy_bits over clones with founder=f)
    S_obs = Var_f(m_f)

Null model:
  - Shuffle founder labels (Starting population) within Time point strata at the cell level.
  - Recompute clone founder assignment as the majority founder within each clone (deterministic).
  - Recompute S_null.

Viability gates (Phase-D):
  - D1: n_clones_total >= min_units_total
  - D2: at least 2 founders have n_clones >= min_units_per_founder
  - D3: perms completed (B >= 1)

Outputs:
  out_dir/
    phaseD_summary.txt
    dispersion_by_founder.tsv
    unit_dispersion_table.tsv
    null_dist_overall.tsv
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _as_boolish(x: object) -> bool:
    if isinstance(x, bool):
        return bool(x)
    if x is None:
        return False
    s = str(x).strip().lower()
    return s in ("1", "true", "t", "yes", "y")


def _entropy_bits_from_counts(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0.0:
        return 0.0
    p = counts.astype(np.float64) / total
    p = p[p > 0.0]
    if p.size == 0:
        return 0.0
    return float(-np.sum(p * np.log2(p)))


def _majority_label(labels: np.ndarray) -> Tuple[object, float]:
    """
    Deterministic majority label + purity.
    purity = max_count / n
    """
    if labels.size == 0:
        return ("", 0.0)
    vc = pd.Series(labels).value_counts(dropna=False)
    top = vc.index[0]
    purity = float(vc.iloc[0] / float(labels.size))
    return (top, purity)


def _p_ge(null: np.ndarray, obs: float) -> float:
    # +1 correction
    B = int(null.size)
    return float((1 + np.sum(null >= obs)) / (1 + B))


@dataclass
class PhaseDConfig:
    founder_col: str = "Starting population"
    time_col: str = "Time point"
    clone_col: str = "clone_id"
    clone_valid_col: str = "clone_valid"


def load_larry_phase1(path_tsv: str, cfg: PhaseDConfig) -> pd.DataFrame:
    df = pd.read_csv(path_tsv, sep="\t")
    required = [cfg.founder_col, cfg.time_col, cfg.clone_col, cfg.clone_valid_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns in input TSV: {missing}")

    # clone_valid is boolean-ish
    df = df.copy()
    df[cfg.clone_valid_col] = df[cfg.clone_valid_col].map(_as_boolish)

    # Keep only valid clones (audit-defined)
    df = df[df[cfg.clone_valid_col] == True].copy()

    # Drop missing clone/time/founder
    df = df[df[cfg.clone_col].notna()].copy()
    df = df[df[cfg.time_col].notna()].copy()
    df = df[df[cfg.founder_col].notna()].copy()

    # Normalize to strings for stable grouping
    df[cfg.clone_col] = df[cfg.clone_col].astype(str)
    df[cfg.time_col] = df[cfg.time_col].astype(str)
    df[cfg.founder_col] = df[cfg.founder_col].astype(str)

    return df


def compute_clone_dispersion_table(df: pd.DataFrame, cfg: PhaseDConfig) -> pd.DataFrame:
    """
    Returns per-clone table:
      clone_id, founder_majority, founder_purity, n_cells, time_entropy_bits, n_time_bins
    """
    rows: List[Dict[str, object]] = []

    for clone_id, g in df.groupby(cfg.clone_col, sort=False):
        times = g[cfg.time_col].to_numpy(dtype=object)
        founders = g[cfg.founder_col].to_numpy(dtype=object)

        founder_major, founder_purity = _majority_label(founders)

        vc_time = pd.Series(times).value_counts()
        ent_bits = _entropy_bits_from_counts(vc_time.to_numpy(dtype=np.int64))

        rows.append(
            dict(
                clone_id=str(clone_id),
                founder_majority=str(founder_major),
                founder_purity=float(founder_purity),
                n_cells=int(g.shape[0]),
                n_time_bins=int(vc_time.shape[0]),
                time_entropy_bits=float(ent_bits),
            )
        )

    out = pd.DataFrame(rows)
    return out


def compute_founder_summary(
    unit_tbl: pd.DataFrame,
    *,
    min_units_per_founder: int,
) -> Tuple[pd.DataFrame, float]:
    """
    Returns:
      founder_summary table with per-founder medians and counts,
      and S = Var_f(median_entropy_bits) over eligible founders.
    """
    rows = []
    eligible_medians = []

    for founder, g in unit_tbl.groupby("founder_majority", sort=True):
        n_units = int(g.shape[0])
        med = float(np.median(g["time_entropy_bits"].to_numpy(dtype=np.float64)))
        rows.append(
            dict(
                founder=str(founder),
                n_clones=n_units,
                median_time_entropy_bits=med,
                mean_time_entropy_bits=float(np.mean(g["time_entropy_bits"].to_numpy(dtype=np.float64))),
                median_founder_purity=float(np.median(g["founder_purity"].to_numpy(dtype=np.float64))),
            )
        )
        if n_units >= int(min_units_per_founder):
            eligible_medians.append(med)

    summ = pd.DataFrame(rows).sort_values(by=["n_clones", "founder"], ascending=[False, True])
    if len(eligible_medians) >= 2:
        S = float(np.var(np.asarray(eligible_medians, dtype=np.float64), ddof=0))
    else:
        S = float("nan")
    return summ, S


def shuffle_founders_within_time(df: pd.DataFrame, cfg: PhaseDConfig, seed: int) -> pd.DataFrame:
    """
    Returns a copy of df with cfg.founder_col shuffled within cfg.time_col strata.
    Clone membership, time labels, and all else preserved.
    """
    rng = np.random.default_rng(seed)
    out = df.copy()

    # Shuffle founder within each time point stratum
    for tp, idx in out.groupby(cfg.time_col, sort=False).groups.items():
        idx = np.asarray(list(idx), dtype=np.int64)
        if idx.size <= 1:
            continue
        out.loc[idx, cfg.founder_col] = rng.permutation(out.loc[idx, cfg.founder_col].to_numpy(dtype=object))

    return out


def write_summary_txt(path: str, lines: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for ln in lines:
            f.write(ln.rstrip() + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_tsv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)

    # Viability gates (fixed defaults; change only with explicit justification)
    ap.add_argument("--min_units_total", type=int, default=200)
    ap.add_argument("--min_units_per_founder", type=int, default=30)

    args = ap.parse_args()
    cfg = PhaseDConfig()

    _ensure_dir(args.out_dir)

    df = load_larry_phase1(args.input_tsv, cfg)

    # Build per-clone unit table (observed)
    unit_tbl = compute_clone_dispersion_table(df, cfg)

    n_units_total = int(unit_tbl.shape[0])
    founders = unit_tbl["founder_majority"].value_counts().to_dict()
    n_founders_total = int(len(founders))
    n_founders_eligible = int(sum(1 for k, v in founders.items() if int(v) >= int(args.min_units_per_founder)))

    # Phase-D viability
    viability_ok = True
    reasons = []

    if n_units_total < int(args.min_units_total):
        viability_ok = False
        reasons.append(f"D1 failed: n_clones_total < {int(args.min_units_total)} (got {n_units_total})")

    if n_founders_eligible < 2:
        viability_ok = False
        reasons.append(
            f"D2 failed: need >=2 founders with n_clones >= {int(args.min_units_per_founder)} "
            f"(got {n_founders_eligible}; founders={n_founders_total})"
        )

    B = int(args.perms)
    if B < 1:
        viability_ok = False
        reasons.append("D3 failed: perms < 1")

    # Observed founder summary + S_obs
    founder_summ_obs, S_obs = compute_founder_summary(
        unit_tbl,
        min_units_per_founder=int(args.min_units_per_founder),
    )

    # Null
    S_null = np.full(B, np.nan, dtype=np.float64)

    if viability_ok:
        for b in range(B):
            df_sh = shuffle_founders_within_time(df, cfg, seed=int(args.seed) + 1000 + b)
            unit_tbl_sh = compute_clone_dispersion_table(df_sh, cfg)
            _, S_b = compute_founder_summary(
                unit_tbl_sh,
                min_units_per_founder=int(args.min_units_per_founder),
            )
            S_null[b] = S_b

        # Some permutations may yield <2 eligible founders (rare); drop NaNs for p-value computation
        S_null_clean = S_null[np.isfinite(S_null)]
        if S_null_clean.size < max(10, int(0.8 * B)):
            # Too many degenerate nulls -> treat as viability failure
            viability_ok = False
            reasons.append(
                f"D3 failed: too many degenerate null perms (finite {S_null_clean.size}/{B})"
            )
            p_S = float("nan")
            S_null_med = float("nan")
            S_null_95 = float("nan")
        else:
            S_null_med = float(np.median(S_null_clean))
            S_null_95 = float(np.quantile(S_null_clean, 0.95))
            p_S = _p_ge(S_null_clean, float(S_obs)) if np.isfinite(S_obs) else float("nan")
    else:
        p_S = float("nan")
        S_null_clean = np.array([], dtype=np.float64)
        S_null_med = float("nan")
        S_null_95 = float("nan")

    # Decision
    if not viability_ok:
        decision = "NO-GO (viability)"
    else:
        # Detectability GO if S_obs > S_null_95
        detect_go = bool(np.isfinite(S_obs) and (float(S_obs) > float(S_null_95)))
        decision = "GO" if detect_go else "NO-GO (detectability)"

    # Write artifacts
    unit_path = os.path.join(args.out_dir, "unit_dispersion_table.tsv")
    unit_tbl.to_csv(unit_path, sep="\t", index=False)

    founder_path = os.path.join(args.out_dir, "dispersion_by_founder.tsv")
    founder_summ_obs.to_csv(founder_path, sep="\t", index=False)

    null_path = os.path.join(args.out_dir, "null_dist_overall.tsv")
    pd.DataFrame({"S_null": S_null}).to_csv(null_path, sep="\t", index=False)

    # Summary
    uplift = (float(S_obs) / float(S_null_med)) if (np.isfinite(S_obs) and np.isfinite(S_null_med) and S_null_med > 0) else float("nan")
    frac_pure = float(np.mean(unit_tbl["founder_purity"].to_numpy(dtype=np.float64) >= 0.999999)) if n_units_total > 0 else float("nan")

    summary_lines = []
    summary_lines.append("=== Phase-D (LARRY) — Founder dispersion / transition admissibility ===")
    summary_lines.append(f"input_tsv={args.input_tsv}")
    summary_lines.append(f"out_dir={args.out_dir}")
    summary_lines.append("")
    summary_lines.append("--- Units (clones) ---")
    summary_lines.append(f"n_clones_total={n_units_total}")
    summary_lines.append(f"n_founders_total={n_founders_total}")
    summary_lines.append(f"n_founders_eligible(min_units_per_founder={int(args.min_units_per_founder)})={n_founders_eligible}")
    summary_lines.append(f"founder_purity_fraction_eq1={frac_pure:.6f}")
    summary_lines.append("")
    summary_lines.append("--- Viability ---")
    summary_lines.append(f"min_units_total={int(args.min_units_total)}")
    summary_lines.append(f"min_units_per_founder={int(args.min_units_per_founder)}")
    summary_lines.append(f"perms={int(args.perms)}")
    summary_lines.append(f"viability_ok={str(viability_ok).lower()}")
    if reasons:
        summary_lines.append("viability_reasons:")
        for r in reasons:
            summary_lines.append(f"  - {r}")
    summary_lines.append("")
    summary_lines.append("--- Observed statistic ---")
    summary_lines.append("metric=Var_f(median_time_entropy_bits_by_founder) over eligible founders")
    summary_lines.append(f"S_obs={S_obs if np.isfinite(S_obs) else 'nan'}")
    summary_lines.append("")
    summary_lines.append("--- Null ---")
    summary_lines.append("null=shuffle Starting population within Time point strata (cell-level), recompute clone majority founder")
    summary_lines.append(f"S_null_median={S_null_med if np.isfinite(S_null_med) else 'nan'}")
    summary_lines.append(f"S_null_95={S_null_95 if np.isfinite(S_null_95) else 'nan'}")
    summary_lines.append(f"uplift_S=S_obs/median(S_null)={uplift if np.isfinite(uplift) else 'nan'}")
    summary_lines.append(f"p_S(P(null>=obs))={p_S if np.isfinite(p_S) else 'nan'}")
    summary_lines.append("")
    summary_lines.append("--- Decision ---")
    summary_lines.append(f"Decision={decision}")
    summary_lines.append("")
    summary_lines.append("Wrote:")
    summary_lines.append(f"  - {os.path.basename(unit_path)}")
    summary_lines.append(f"  - {os.path.basename(founder_path)}")
    summary_lines.append(f"  - {os.path.basename(null_path)}")

    summary_path = os.path.join(args.out_dir, "phaseD_summary.txt")
    write_summary_txt(summary_path, summary_lines)

    # Console echo (concise)
    print("=== Phase-D summary ===")
    print(f"n_clones_total={n_units_total}")
    print(f"n_founders_eligible={n_founders_eligible}")
    print(f"viability_ok={str(viability_ok).lower()}")
    if reasons:
        for r in reasons:
            print(f"  {r}")
    print(f"S_obs={S_obs if np.isfinite(S_obs) else 'nan'}")
    print(f"S_null_median={S_null_med if np.isfinite(S_null_med) else 'nan'}")
    print(f"S_null_95={S_null_95 if np.isfinite(S_null_95) else 'nan'}")
    print(f"p_S={p_S if np.isfinite(p_S) else 'nan'}")
    print(f"Decision={decision}")
    print(f"Wrote: {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

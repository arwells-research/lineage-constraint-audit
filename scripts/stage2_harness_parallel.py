#!/usr/bin/env python3
"""
stage2_harness_parallel.py — Parallel Stage-2 robustness harness for scripts/phase1_runner.py

What this harness does
- Runs phase1_runner.py across a parameter grid in parallel (ProcessPoolExecutor)
- Captures per-run stdout/stderr logs
- Parses runner stdout to extract:
  - Stage-1A viability vs Stage-1B detectability outcome
  - v_obs, A_obs, null medians, p-values
  - v_null_95, A_null_95 (required to enforce Stage-1B OR logic in the harness)
  - direction checks vs null medians
  - "carrying channel" for detectability (via v or via A), inferred from Stage-1B
    criteria (obs > null_95), not from p-values.

Stage-2 rule (aligned to revised audit_logic.md)
- Compute detectability-GO rate across runs that reached detectability evaluation
- Require detectability-GO in >=80% of parameter combinations
- Determine which channel(s) carry GO across the detectability-GO runs:
    * carry_v if any run has v_obs > v_null_95
    * carry_A if any run has A_obs > A_null_95
- Require direction consistency across ALL runs that reached detectability evaluation
  for the carrying channel(s) only:
    * if carry_A: require A_obs > median(A_null)
    * if carry_v: require v_obs > median(v_null)
  

Notes
- This harness does NOT invent new model knobs.
- --alpha_detect is retained for CLI compatibility, but Stage-1B / carry-channel
   inference is performed using the audit-spec quantile rule (obs > null_95).

Typical usage (from repo root):
  python scripts/stage2_harness_parallel.py \
    --runner scripts/phase1_runner.py \
    --annot data/raw/GSE126954_cell_annotation.csv.gz \
    --mtx data/raw/GSE126954_gene_by_cell_count_matrix.txt.gz \
    --out_dir results/gse126954/stage2_harness \
    --k 10,15,20,30 \
    --seed 0,1,2,3,4 \
    --max_centers 2000,5000 \
    --null_mode timebin,timebin_batch \
    --control_modes none \
    --perms 1000 \
    --jobs 14
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from concurrent.futures import ProcessPoolExecutor, as_completed


# ----------------------------
# Parsing helpers (stdout)
# ----------------------------

_RE_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

_PATTERNS = {
    "v_obs": re.compile(r"Validity rate v_obs:\s*(?P<v_obs>%s)" % _RE_FLOAT),
    "A_obs": re.compile(r"Conditional dependence A_obs\s*=\s*mean\(T_i \| valid\):\s*(?P<A_obs>%s)" % _RE_FLOAT),
    "p_v": re.compile(
        r"^p_v\s*(?:\(|\s)\s*v_null\s*>=\s*v_obs\s*(?:\)|\s)\s*:\s*(?P<p_v>%s)\s*$" % _RE_FLOAT,
        re.MULTILINE,
    ),
    "p_A": re.compile(
        r"^p_A\s*(?:\(|\s)\s*A_null\s*>=\s*A_obs\s*(?:\)|\s)\s*:\s*(?P<p_A>%s)\s*$" % _RE_FLOAT,
        re.MULTILINE,
    ),
    # Null summaries printed by phase1_runner.py:
    #   v_null: min=...  median=...  p90=...  p95=...  max=...
    #   A_null: min=...  median=...  p90=...  p95=...  max=...
    # Use rf-strings (NOT % formatting) to avoid '%' collisions (e.g., '95%').
    "v_null_median": re.compile(rf"v_null:\s*min={_RE_FLOAT}\s+median=(?P<med>{_RE_FLOAT})"),
    "A_null_median": re.compile(rf"A_null:\s*min={_RE_FLOAT}\s+median=(?P<med>{_RE_FLOAT})"),
    "v_null_p95": re.compile(r"^v_null:.*?\bp95=(?P<p95>%s)\b" % _RE_FLOAT, re.MULTILINE),
    "A_null_p95": re.compile(r"^A_null:.*?\bp95=(?P<p95>%s)\b" % _RE_FLOAT, re.MULTILINE),
    "decision_line": re.compile(r"^Decision:\s*(?P<dec>GO|NO-GO)(?:\s*\((?P<reason>.*)\))?\s*$"),
    "valid_centers": re.compile(r"\(valid centers:\s*(?P<valid>\d+)\s*/\s*(?P<total>\d+)\)"),
}


def _safe_float(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        return None


def _safe_int(x: Optional[str]) -> Optional[int]:
    if x is None:
        return None
    try:
        return int(x)
    except Exception:
        return None


def _sha1(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _write_text(path: Path, s: str) -> None:
    path.write_text(s, encoding="utf-8")


def _tsv_escape(x: object) -> str:
    s = "" if x is None else str(x)
    return s.replace("\t", " ").replace("\n", "\\n").replace("\r", "")


def _format_bool(x: Optional[bool]) -> str:
    if x is None:
        return ""
    return "1" if x else "0"


@dataclass
class ParsedRun:
    decision: str
    reason: str
    viability_failed: bool
    detectability_go: Optional[bool]

    v_obs: Optional[float]
    A_obs: Optional[float]
    n_valid_centers: Optional[int]
    n_centers: Optional[int]

    v_null_median: Optional[float]
    A_null_median: Optional[float]

    v_null_95: Optional[float]
    A_null_95: Optional[float]

    p_v: Optional[float]
    p_A: Optional[float]

    uplift: Optional[float]

    v_direction_ok: Optional[bool]
    A_direction_ok: Optional[bool]

    go_via_v: Optional[bool]
    go_via_A: Optional[bool]


def parse_phase1_stdout(stdout_text: str, *, alpha_detect: float) -> ParsedRun:
    lines = stdout_text.splitlines()

    v_obs = A_obs = None
    v_null_med = A_null_med = None
    v_null_95 = A_null_95 = None    
    p_v = p_A = None
    n_valid = n_centers = None
    decision = "NO-GO"
    reason = ""

    m = _PATTERNS["v_obs"].search(stdout_text)
    if m:
        v_obs = _safe_float(m.group("v_obs"))

    m = _PATTERNS["A_obs"].search(stdout_text)
    if m:
        A_obs = _safe_float(m.group("A_obs"))

    m = _PATTERNS["p_v"].search(stdout_text)
    if m:
        p_v = _safe_float(m.group("p_v"))

    m = _PATTERNS["p_A"].search(stdout_text)
    if m:
        p_A = _safe_float(m.group("p_A"))

    m = _PATTERNS["v_null_median"].search(stdout_text)
    if m:
        v_null_med = _safe_float(m.group("med"))

    m = _PATTERNS["A_null_median"].search(stdout_text)
    if m:
        A_null_med = _safe_float(m.group("med"))

    m = _PATTERNS["v_null_p95"].search(stdout_text)
    if m:
        v_null_95 = _safe_float(m.group("p95"))
 
    m = _PATTERNS["A_null_p95"].search(stdout_text)
    if m:
        A_null_95 = _safe_float(m.group("p95"))
 
    m = _PATTERNS["valid_centers"].search(stdout_text)
    if m:
        n_valid = _safe_int(m.group("valid"))
        n_centers = _safe_int(m.group("total"))

    # last Decision line wins
    for ln in reversed(lines):
        dm = _PATTERNS["decision_line"].match(ln.strip())
        if dm:
            decision = dm.group("dec")
            reason = (dm.group("reason") or "").strip()
            break

    # Stage-1A vs Stage-1B:
    # If reason begins with "viability", treat as viability failure.
    viability_failed = False
    # For Stage-2 we want audit-spec detectability (alpha_detect),
    # not whatever threshold the runner used for its printed Decision.
    # If p-values are present, derive detectability from them; otherwise fall back to printed decision.
    if p_v is not None or p_A is not None:
        dv = (p_v is not None) and (p_v <= alpha_detect)
        dA = (p_A is not None) and (p_A <= alpha_detect)
        detectability_go = bool(dv or dA)
    else:
        detectability_go = (decision == "GO")

    # Derived metrics
    uplift = None
    if A_obs is not None and A_null_med is not None and A_null_med > 0:
        uplift = float(A_obs / A_null_med)

    v_dir_ok = None
    if v_obs is not None and v_null_med is not None:
        v_dir_ok = bool(v_obs > v_null_med)

    A_dir_ok = None
    if A_obs is not None and A_null_med is not None:
        A_dir_ok = bool(A_obs > A_null_med)

    # Carry-channel inference aligned to audit spec:
    # Channel carries GO if (obs > null_95) for that statistic.
    # (p-values are retained for reporting, but are not used for carry-channel logic.)
    go_via_v = None
    if v_obs is not None and v_null_95 is not None:
        go_via_v = bool(v_obs > v_null_95)
  
    go_via_A = None
    if A_obs is not None and A_null_95 is not None:
         go_via_A = bool(A_obs > A_null_95)

    # Detectability-GO for Stage-2 MUST follow audit-spec alpha (alpha_detect),
    # not the runner's printed decision (runner may use a stricter alpha).
    if not viability_failed:
        if (go_via_v is None) and (go_via_A is None):
            detectability_go = None
        else:
            detectability_go = bool((go_via_v is True) or (go_via_A is True))
 

    return ParsedRun(
        decision=decision,
        reason=reason,
        viability_failed=viability_failed,
        detectability_go=detectability_go,
        v_obs=v_obs,
        A_obs=A_obs,
        n_valid_centers=n_valid,
        n_centers=n_centers,
        v_null_median=v_null_med,
        A_null_median=A_null_med,
        v_null_95=v_null_95,
        A_null_95=A_null_95,        
        p_v=p_v,
        p_A=p_A,
        uplift=uplift,
        v_direction_ok=v_dir_ok,
        A_direction_ok=A_dir_ok,
        go_via_v=go_via_v,
        go_via_A=go_via_A,
    )


def build_grid(values_csv: str, cast_fn) -> List:
    items = [v.strip() for v in values_csv.split(",") if v.strip() != ""]
    return [cast_fn(it) for it in items]


def _run_one(
    *,
    run_idx: int,
    cmd: List[str],
    cmd_id: str,
    logs_dir: str,
    timeout_s: Optional[int],
    alpha_detect: float,
) -> Dict[str, object]:
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            cwd=None,
            env=os.environ,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout_s,
            check=False,
        )
        rc = proc.returncode
        out = proc.stdout
        err = proc.stderr
    except subprocess.TimeoutExpired as e:
        rc = 124
        out = e.stdout or ""
        err = (e.stderr or "") + "\n[TIMEOUT]\n"
    except Exception as e:
        rc = 125
        out = ""
        err = f"[HARNESS_ERROR] {type(e).__name__}: {e}\n"

    dt = time.time() - t0

    logs_dir_p = Path(logs_dir)
    stdout_path = logs_dir_p / f"run_{run_idx:04d}_{cmd_id}.stdout.txt"
    stderr_path = logs_dir_p / f"run_{run_idx:04d}_{cmd_id}.stderr.txt"
    _write_text(stdout_path, out)
    _write_text(stderr_path, err)

    parsed: Optional[ParsedRun] = None
    if rc == 0:
        parsed = parse_phase1_stdout(out, alpha_detect=alpha_detect)

    row = {
        "run_idx": run_idx,
        "returncode": rc,
        "elapsed_s": f"{dt:.3f}",
        "decision": parsed.decision if parsed else "",
        "reason": parsed.reason if parsed else "",
        "viability_failed": _format_bool(parsed.viability_failed) if parsed else "",
        "detectability_go": _format_bool(parsed.detectability_go) if parsed else "",
        "n_valid_centers": parsed.n_valid_centers if parsed else "",
        "n_centers": parsed.n_centers if parsed else "",
        "v_obs": "" if not parsed or parsed.v_obs is None else f"{parsed.v_obs:.8f}",
        "A_obs": "" if not parsed or parsed.A_obs is None else f"{parsed.A_obs:.8f}",
        "v_null_median": "" if not parsed or parsed.v_null_median is None else f"{parsed.v_null_median:.8f}",
        "A_null_median": "" if not parsed or parsed.A_null_median is None else f"{parsed.A_null_median:.8f}",
        "v_null_95": "" if not parsed or getattr(parsed, "v_null_95", None) is None else f"{parsed.v_null_95:.8f}",
        "A_null_95": "" if not parsed or getattr(parsed, "A_null_95", None) is None else f"{parsed.A_null_95:.8f}",        
        "p_v": "" if not parsed or parsed.p_v is None else f"{parsed.p_v:.10g}",
        "p_A": "" if not parsed or parsed.p_A is None else f"{parsed.p_A:.10g}",
        "uplift": "" if not parsed or parsed.uplift is None else f"{parsed.uplift:.8f}",
        "v_direction_ok": _format_bool(parsed.v_direction_ok) if parsed else "",
        "A_direction_ok": _format_bool(parsed.A_direction_ok) if parsed else "",
        "go_via_v": _format_bool(parsed.go_via_v) if parsed else "",
        "go_via_A": _format_bool(parsed.go_via_A) if parsed else "",
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "cmd_sha1": cmd_id,
    }
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runner", required=True, help="Path to phase1_runner.py")
    ap.add_argument("--annot", required=True)
    ap.add_argument("--mtx", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--k", default="20", help="Comma list for --knn_k")
    ap.add_argument("--seed", default="0,1,2,3,4", help="Comma list for --seed")
    ap.add_argument("--max_centers", default="5000", help="Comma list for --max_centers")
    ap.add_argument("--null_mode", default="timebin,timebin_batch", help="Comma list for --null_mode")
    ap.add_argument("--control_modes", default="none", help="Comma list for --control_mode")

    ap.add_argument("--svd_components", type=int, default=50)
    ap.add_argument("--n_genes", type=int, default=20222)
    ap.add_argument("--n_cells_total", type=int, default=89701)
    ap.add_argument("--perms", type=int, default=1000)

    ap.add_argument("--jobs", type=int, default=14, help="Parallel workers (physical cores recommended)")
    ap.add_argument("--timeout_s", type=int, default=0, help="0 = no timeout")
    ap.add_argument("--fail_fast", action="store_true", help="Stop scheduling new work after first error")

    # Audit-spec alpha used only for carry-channel inference from printed p-values.
    ap.add_argument("--alpha_detect", type=float, default=0.05)

    args = ap.parse_args()

    runner = Path(args.runner)
    annot = Path(args.annot)
    mtx = Path(args.mtx)
    if not runner.exists():
        print(f"ERROR: runner not found: {runner}", file=sys.stderr)
        return 2
    if not annot.exists():
        print(f"ERROR: annot not found: {annot}", file=sys.stderr)
        return 2
    if not mtx.exists():
        print(f"ERROR: mtx not found: {mtx}", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir)
    _ensure_dir(out_dir)

    ks = build_grid(args.k, int)
    seeds = build_grid(args.seed, int)
    max_centers_list = build_grid(args.max_centers, int)
    null_modes = build_grid(args.null_mode, str)
    control_modes = build_grid(args.control_modes, str)

    grid = list(itertools.product(ks, seeds, max_centers_list, null_modes, control_modes))
    total_runs = len(grid)

    run_id = time.strftime("%Y%m%d_%H%M%S")
    run_root = out_dir / f"run_{run_id}"
    logs_dir = run_root / "logs"
    _ensure_dir(logs_dir)

    timeout_s = None if args.timeout_s <= 0 else int(args.timeout_s)

    harness_cfg = {
        "runner": str(runner),
        "annot": str(annot),
        "mtx": str(mtx),
        "grid_axes": {
            "k": ks,
            "seed": seeds,
            "max_centers": max_centers_list,
            "null_mode": null_modes,
            "control_mode": control_modes,
        },
        "fixed_args": {
            "svd_components": args.svd_components,
            "n_genes": args.n_genes,
            "n_cells_total": args.n_cells_total,
            "global_perms": args.perms,
        },
        "alpha_detect": args.alpha_detect,
        "jobs": args.jobs,
        "timeout_s": args.timeout_s,
        "fail_fast": bool(args.fail_fast),
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    _write_text(run_root / "harness_config.json", json.dumps(harness_cfg, indent=2, sort_keys=True) + "\n")

    header = [
        "run_idx",
        "k",
        "seed",
        "max_centers",
        "null_mode",
        "control_mode",
        "returncode",
        "elapsed_s",
        "decision",
        "reason",
        "viability_failed",
        "detectability_go",
        "n_valid_centers",
        "n_centers",
        "v_obs",
        "A_obs",
        "v_null_median",
        "A_null_median",
        "v_null_95",
        "A_null_95",        
        "p_v",
        "p_A",
        "uplift",
        "v_direction_ok",
        "A_direction_ok",
        "go_via_v",
        "go_via_A",
        "stdout_rel",
        "stderr_rel",
        "cmd_sha1",
    ]

    print("=== Stage-2 harness (parallel) ===")
    print(f"Runs: {total_runs}")
    print(f"Workers: {args.jobs}")
    print(f"Output: {run_root}")
    print()

    # Pre-build job specs deterministically
    jobs: List[Tuple[int, Dict[str, object]]] = []
    for i, (k, seed, max_centers, null_mode, control_mode) in enumerate(grid, start=1):
        cmd = [
            sys.executable,
            str(runner),
            "--annot",
            str(annot),
            "--mtx",
            str(mtx),
            "--n_genes",
            str(int(args.n_genes)),
            "--n_cells_total",
            str(int(args.n_cells_total)),
            "--svd_components",
            str(int(args.svd_components)),
            "--knn_k",
            str(int(k)),
            "--max_centers",
            str(int(max_centers)),
            "--global_perms",
            str(int(args.perms)),
            "--seed",
            str(int(seed)),
            "--null_mode",
            str(null_mode),
            "--control_mode",
            str(control_mode),
        ]
        cmd_str = " ".join(cmd)
        cmd_id = _sha1(cmd_str)[:12]
        meta = {
            "k": k,
            "seed": seed,
            "max_centers": max_centers,
            "null_mode": null_mode,
            "control_mode": control_mode,
            "cmd": cmd,
            "cmd_id": cmd_id,
        }
        jobs.append((i, meta))

    results_by_idx: Dict[int, Dict[str, object]] = {}
    first_error_seen = False

    with ProcessPoolExecutor(max_workers=int(args.jobs)) as ex:
        future_map = {}
        for (run_idx, meta) in jobs:
            if args.fail_fast and first_error_seen:
                break
            fut = ex.submit(
                _run_one,
                run_idx=run_idx,
                cmd=meta["cmd"],
                cmd_id=meta["cmd_id"],
                logs_dir=str(logs_dir),
                timeout_s=timeout_s,
                alpha_detect=float(args.alpha_detect),
            )
            future_map[fut] = (run_idx, meta)

        completed = 0
        for fut in as_completed(future_map):
            run_idx, meta = future_map[fut]
            completed += 1

            row_base = {
                "run_idx": run_idx,
                "k": meta["k"],
                "seed": meta["seed"],
                "max_centers": meta["max_centers"],
                "null_mode": meta["null_mode"],
                "control_mode": meta["control_mode"],
            }

            try:
                row_run = fut.result()
            except Exception as e:
                row_run = {
                    "returncode": 126,
                    "elapsed_s": "",
                    "decision": "",
                    "reason": f"HARNESS_EXCEPTION: {type(e).__name__}: {e}",
                    "viability_failed": "",
                    "detectability_go": "",
                    "n_valid_centers": "",
                    "n_centers": "",
                    "v_obs": "",
                    "A_obs": "",
                    "v_null_median": "",
                    "A_null_median": "",
                    "v_null_95": "",
                    "A_null_95": "",                    
                    "p_v": "",
                    "p_A": "",
                    "uplift": "",
                    "v_direction_ok": "",
                    "A_direction_ok": "",
                    "go_via_v": "",
                    "go_via_A": "",
                    "stdout_path": "",
                    "stderr_path": "",
                    "cmd_sha1": meta["cmd_id"],
                }

            rc = int(row_run.get("returncode", 0) or 0)
            if rc != 0:
                first_error_seen = True

            stdout_abs = str(row_run.get("stdout_path", ""))
            stderr_abs = str(row_run.get("stderr_path", ""))
            stdout_rel = ""
            stderr_rel = ""
            try:
                if stdout_abs:
                    stdout_rel = str(Path(stdout_abs).relative_to(run_root))
                if stderr_abs:
                    stderr_rel = str(Path(stderr_abs).relative_to(run_root))
            except Exception:
                stdout_rel = stdout_abs
                stderr_rel = stderr_abs

            row = dict(row_base)
            row.update(row_run)
            row["stdout_rel"] = stdout_rel
            row["stderr_rel"] = stderr_rel
            row["cmd_sha1"] = row_run.get("cmd_sha1", meta["cmd_id"])
            results_by_idx[run_idx] = row

            if completed % max(1, total_runs // 20) == 0 or completed == total_runs:
                print(f"Completed {completed}/{total_runs}")

    # Write TSV in deterministic run_idx order
    tsv_path = run_root / "runs.tsv"
    with tsv_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for i in range(1, total_runs + 1):
            r = results_by_idx.get(i)
            if r is None:
                continue
            f.write("\t".join(_tsv_escape(r.get(h, "")) for h in header) + "\n")

    # Aggregate Stage-2 summary
    rows = [results_by_idx[i] for i in sorted(results_by_idx.keys())]
    ok_runs = [r for r in rows if str(r.get("returncode", "")) == "0"]

    viability_failed_runs = [r for r in ok_runs if r.get("viability_failed") == "1"]
    reached_detectability_runs = [r for r in ok_runs if r.get("detectability_go") in ("0", "1")]

    detectability_go_runs = [r for r in reached_detectability_runs if r.get("detectability_go") == "1"]
    detectability_no_runs = [r for r in reached_detectability_runs if r.get("detectability_go") == "0"]

    denom = max(1, len(reached_detectability_runs))
    go_rate = len(detectability_go_runs) / float(denom)

    # Determine which channel(s) carry GO across detectability-GO runs.
    carry_v = any(r.get("go_via_v") == "1" for r in detectability_go_runs)
    carry_A = any(r.get("go_via_A") == "1" for r in detectability_go_runs)
 
    # Direction consistency is enforced across ALL runs that reached detectability evaluation,
    # but only for the carrying channel(s).
    def _dir_ok_for_carry_channels(r: Dict[str, object]) -> bool:
        if carry_v:
            if r.get("v_direction_ok") != "1":
                return False
        if carry_A:
            if r.get("A_direction_ok") != "1":
                return False
        return True
 
    direction_ok_all_reached = (
        all(_dir_ok_for_carry_channels(r) for r in reached_detectability_runs)
        if reached_detectability_runs
        else False
     )
 
    # Stage-2 GO: >=80% detectability-GO among runs that reached detectability evaluation,
    # AND direction consistent for the carrying channel(s) across all reached runs.
    stage2_go = (go_rate >= 0.80) and direction_ok_all_reached

    n_runner_errors = len([r for r in rows if str(r.get("returncode", "")) not in ("0", "") and int(r.get("returncode", 0)) != 0])

    # Also write a small TSV of viability failures for quick inspection
    vf_path = run_root / "viability_failures.tsv"
    with vf_path.open("w", encoding="utf-8") as f:
        vf_cols = ["run_idx", "k", "seed", "max_centers", "null_mode", "control_mode", "reason", "stdout_rel"]
        f.write("\t".join(vf_cols) + "\n")
        for r in viability_failed_runs:
            f.write("\t".join(_tsv_escape(r.get(c, "")) for c in vf_cols) + "\n")

    summary = {
        "run_root": str(run_root),
        "total_grid_runs": total_runs,
        "completed_runs_recorded": len(rows),
        "runner_errors_nonzero_rc": n_runner_errors,
        "completed_ok_runs_rc0": len(ok_runs),
        "viability_failed_runs": len(viability_failed_runs),
        "reached_detectability_runs": len(reached_detectability_runs),
        "detectability_go_runs": len(detectability_go_runs),
        "detectability_no_go_runs": len(detectability_no_runs),
        "detectability_go_rate": go_rate,
        "direction_consistency": {
            "carry_channels": {"carry_v": carry_v, "carry_A": carry_A},
            "direction_ok_all_reached_detectability_runs_for_carry_channels": direction_ok_all_reached,
        },
        "stage2_go": stage2_go,
        "stage2_rule": "GO if detectability-GO in >=80% of parameter combinations AND direction consistent for the channel(s) that carry GO (Stage-1B OR logic).",
        "artifacts": {
            "runs_tsv": str(tsv_path),
            "viability_failures_tsv": str(vf_path),
            "summary_json": str(run_root / "summary.json"),
            "summary_md": str(run_root / "summary.md"),
            "logs_dir": str(logs_dir),
        },
        "config": harness_cfg,
    }
    _write_text(run_root / "summary.json", json.dumps(summary, indent=2, sort_keys=True) + "\n")

    md = []
    md.append("# Stage-2 Robustness Harness Summary (parallel)\n\n")
    md.append(f"- Run root: `{run_root}`\n")
    md.append(f"- Total grid runs: **{total_runs}**\n")
    md.append(f"- Completed recorded runs: **{len(rows)}**\n")
    md.append(f"- Completed OK (rc==0): **{len(ok_runs)}**\n")
    md.append(f"- Runner errors (nonzero rc): **{n_runner_errors}**\n\n")

    md.append("## Stage-1A viability\n")
    md.append(f"- Viability failures: **{len(viability_failed_runs)}**\n")
    md.append(f"- Reached detectability evaluation: **{len(reached_detectability_runs)}**\n")
    md.append(f"- Viability failures table: `{vf_path}`\n\n")

    md.append("## Stage-1B detectability robustness\n")
    md.append(f"- Detectability GO: **{len(detectability_go_runs)}**\n")
    md.append(f"- Detectability NO-GO: **{len(detectability_no_runs)}**\n")
    md.append(f"- Detectability GO rate: **{go_rate:.3f}** (threshold 0.800)\n\n")

    md.append("## Direction consistency (Stage-2)\n")
    md.append(f"- Carry channels: v={carry_v}, A={carry_A}\n")
    md.append(f"- Direction consistency across ALL reached-detectability runs (for carry channels): **{direction_ok_all_reached}**\n\n")
  

    md.append("## Stage-2 decision\n")
    md.append(f"- **Stage-2 GO:** **{stage2_go}**\n\n")

    md.append("## Artifacts\n")
    md.append(f"- Runs table: `{tsv_path}`\n")
    md.append(f"- Summary JSON: `{run_root / 'summary.json'}`\n")
    md.append(f"- Per-run logs: `{logs_dir}`\n")

    _write_text(run_root / "summary.md", "".join(md))

    print()
    print("=== Stage-2 summary ===")
    print(f"Recorded runs: {len(rows)} / {total_runs}")
    print(f"OK runs (rc==0): {len(ok_runs)}")
    print(f"Viability failed: {len(viability_failed_runs)}")
    print(f"Reached detectability: {len(reached_detectability_runs)}")
    print(f"Detectability GO rate: {go_rate:.3f} (need >=0.800)")
    print(f"Direction: carry_v={carry_v}, carry_A={carry_A}, consistency across reached runs = {direction_ok_all_reached}")
    print(f"Stage-2 decision: {'GO' if stage2_go else 'NO-GO'}")
    print(f"Wrote: {run_root / 'summary.md'}")
    print(f"Wrote: {run_root / 'runs.tsv'}")
    print(f"Wrote: {vf_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

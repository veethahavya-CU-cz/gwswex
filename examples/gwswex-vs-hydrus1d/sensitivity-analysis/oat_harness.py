# type: ignore
"""Iterated OAT (coordinate-descent) parameter optimisation for GWSWEX.

Subprocess-isolated trials: each candidate parameter set is evaluated by
spawning a fresh `oat_worker.py` process. This protects the harness from
Fortran-kernel state corruption when an individual trial crashes mid-init
("Attempting to allocate already allocated variable 'model'").

For each (soil, setup, solver):

  1. Build a namespace by executing the first 4 code cells of the notebook
     (cell 0 imports, cell 1 config, cell 2 HYDRUS-1D reference, cell 3
     model builder). Cache the HYDRUS WT reference (`h_gw_d`) to
     `oat_results/cache/h_gw_d-{soil}-{setup}.npy` so the worker can skip
     the slow HYDRUS rerun.
  2. Parse MODEL_PARAMS, ET_STRESS, set_solver(...) kwargs from the
     notebook source -> baseline (lets hand-tuned values like basic-clay
     ICratio_min=0.42 survive into OAT).
  3. Coordinate-descent OAT: per pass, sweep one parameter at a time
     (others held at the current best). Accept the per-parameter winner
     if it improves the *phase-targeted metric* by >= IMPROVE_TOL
     relative to the current best on that metric. Repeat up to MAX_PASSES
     or until overall RMSE does not improve by >= CONVERGE_TOL.
  4. Acceptance margin per setup:
        basic     : >= 2 %  (IMPROVE_TOL_BASIC)
        intensive : >= 5 %  (IMPROVE_TOL_INT)

Phase-targeted metrics (PARAM_METRIC table):
  Each parameter is scored against the metric that best reveals its
  physical role. Acceptance of a candidate is based on that metric;
  convergence termination uses overall RMSE.
  - Implicit: picard_tol, picard_max_iter, n_trapz, beta_hyst -> overall
  - Implicit ET stress (s_star, s_w, s_h, s_e)               -> dry
  - Explicit: courant_number, n_trapz                         -> overall
  - Explicit: ICratio_min, F_min, psi_f                       -> wet
  - Explicit: beta_hyst                                       -> dry_cool
  - Explicit ET stress (s_star, s_w, s_h, s_e)               -> dry
  - psi_f (common to both)                                    -> wet

Output: oat_results/oat_results.json
"""

from __future__ import annotations
import argparse
import json
import os
import re
import subprocess
import sys
import time
from datetime import timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.show = lambda: None  # type: ignore

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
EX_DIR = ROOT / "examples" / "gwswex-vs-hydrus1d"
OUT_DIR = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results"
CACHE_DIR = OUT_DIR / "cache"
WORKER = Path(__file__).parent / "oat_worker.py"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

ALL_CASES = ["loam", "sand", "clay", "sand-loam", "sand-clay", "loam-clay"]
ALL_SETUPS = ["basic", "intensive"]

# --- Acceptance margins (relative WT-RMSE improvement) ----------------------
IMPROVE_TOL_BASIC = 0.02
IMPROVE_TOL_INT = 0.05
CONVERGE_TOL = 0.005
MAX_PASSES = 3

WORKER_TIMEOUT_S = 240

# --- Sweep grids ------------------------------------------------------------

EXPLICIT_GRID = {
    "courant_number": [0.3, 0.5, 0.7, 0.85, 0.9, 0.95],
    "n_trapz": [5, 10, 15, 20, 30, 40],
    "beta_hyst": [0.7, 0.85, 0.9, 1.0],
}
IMPLICIT_GRID = {
    "picard_tol": [1e-7, 1e-6, 1e-5, 1e-4],
    "beta_hyst": [0.7, 0.85, 0.9, 1.0],
    "n_trapz": [10, 20, 30, 40],
}
MODEL_GRID = {
    "psi_f": [0.005, 0.01, 0.05, 0.09, 0.15, 0.20],
    "F_min": [1e-8, 1e-7, 1e-6, 1e-5],
    "ICratio_min": [0.05, 0.10, 0.20, 0.30, 0.42, 0.50, 0.60],
}
ET_STRESS_GRID = {
    "s_star": [0.3, 0.4, 0.5, 0.6, 0.7],
    "s_w": [0.05, 0.10, 0.15, 0.20],
    "s_h": [0.02, 0.05, 0.08],
    "s_e": [0.20, 0.30, 0.40, 0.50, 0.60],
}

SOLVER_KEYS = {"courant_number", "n_trapz", "beta_hyst", "picard_tol", "picard_max_iter"}

# ---------------------------------------------------------------------------
# Phase-targeted metric routing
# ---------------------------------------------------------------------------
# Maps (solver, param_name) -> metric_key.
# metric_key is one of: "overall", "wet", "dry", "dry_cool".
# For basic setup, dry_cool == dry (no cooldown phase).
# Parameters not present in either dict default to "overall".
#
# Physical rationale (condensed):
#   "wet"  params: affect infiltration / wetting-front speed; scored where
#          precipitation is active (wet phase).
#   "dry"  params: affect ET partitioning / capillary recession; scored in
#          the dry phase where those processes dominate GWH behaviour.
#   "dry_cool": hysteresis blend; affects the drying branch of the
#          retention curve which spans dry + cooldown in intensive setups.
#   "overall": Picard convergence and quadrature params affect both phases
#          roughly equally; overall RMSE is the appropriate aggregator.
_PARAM_METRIC_IMPLICIT: dict[str, str] = {
    "picard_tol": "overall",
    "picard_max_iter": "overall",
    "n_trapz": "overall",
    "beta_hyst": "overall",
    "s_star": "dry",
    "s_w": "dry",
    "s_h": "dry",
    "s_e": "dry",
}
_PARAM_METRIC_EXPLICIT: dict[str, str] = {
    "courant_number": "overall",
    "n_trapz": "overall",
    "ICratio_min": "wet",
    "F_min": "wet",
    "beta_hyst": "dry_cool",
    "s_star": "dry",
    "s_w": "dry",
    "s_h": "dry",
    "s_e": "dry",
}
_PARAM_METRIC_COMMON: dict[str, str] = {
    "psi_f": "wet",  # Green-Ampt suction head: active only when P > K_sat
}


def _param_metric(solver: str, pname: str) -> str:
    """Return the phase-targeted metric key for this (solver, parameter) pair."""
    if solver == "implicit":
        return _PARAM_METRIC_IMPLICIT.get(pname, _PARAM_METRIC_COMMON.get(pname, "overall"))
    return _PARAM_METRIC_EXPLICIT.get(pname, _PARAM_METRIC_COMMON.get(pname, "overall"))


_NAN_METRICS: dict[str, float] = {
    "overall": float("nan"),
    "wet": float("nan"),
    "dry": float("nan"),
    "dry_cool": float("nan"),
}


def _nb_name(soil_tag: str, setup: str) -> str:
    if setup == "basic":
        nb = EX_DIR / f"comparison-basic-{soil_tag}.ipynb"
        if not nb.exists():
            nb = EX_DIR / f"comparison-{soil_tag}.ipynb"
        return nb.name
    return f"comparison-intensive-{soil_tag}.ipynb"


def _load_cells(nb_path: Path) -> list[dict]:
    return json.loads(nb_path.read_text())["cells"]


def _src(cell) -> str:
    s = cell["source"]
    return "".join(s) if isinstance(s, list) else s


_DICT_KW_RE = re.compile(r"([A-Za-z_]\w*)\s*=\s*([^,]+)")
_SET_SOLVER_RE = re.compile(r'm\.set_solver\(\s*solver\s*=\s*"(?P<solver>explicit|implicit)"(?P<kw>[^)]*)\)')


def _parse_kwargs(blob: str) -> dict:
    out: dict = {}
    for m in _DICT_KW_RE.finditer(blob):
        k = m.group(1)
        v = m.group(2).strip().rstrip(",").strip()
        try:
            out[k] = eval(v, {"__builtins__": {}})
        except Exception:
            out[k] = v
    return out


def _parse_baselines(soil_tag: str, setup: str) -> dict:
    nb = EX_DIR / _nb_name(soil_tag, setup)
    text = "\n".join(_src(c) for c in _load_cells(nb) if c.get("cell_type") == "code")
    out = {"model_params": {}, "et_stress": {}, "explicit_solver_params": {}, "implicit_solver_params": {}}
    m = re.search(r"MODEL_PARAMS\s*=\s*dict\(([^)]*)\)", text)
    if m:
        out["model_params"] = _parse_kwargs(m.group(1))
    m = re.search(r"ET_STRESS\s*=\s*dict\(([^)]*)\)", text)
    if m:
        out["et_stress"] = _parse_kwargs(m.group(1))
    for sm in _SET_SOLVER_RE.finditer(text):
        kw = _parse_kwargs(sm.group("kw").lstrip(",").strip())
        out[f'{sm.group("solver")}_solver_params'] = kw
    return out


def _ensure_hydrus_cache(soil_tag: str, setup: str) -> Path:
    """Run notebook cells 0-2 once, persist h_gw_d to cache, return path."""
    cache = CACHE_DIR / f"h_gw_d-{soil_tag}-{setup}.npy"
    if cache.exists():
        return cache
    nb = EX_DIR / _nb_name(soil_tag, setup)
    cells = _load_cells(nb)
    code_cells = [c for c in cells if c.get("cell_type") == "code"]
    ns: dict = {"__name__": "__oat_cache__", "__file__": str(nb), "__setup__": setup}
    cwd = os.getcwd()
    os.chdir(EX_DIR)
    try:
        for i, c in enumerate(code_cells[:3]):  # 0, 1, 2 (HYDRUS)
            exec(compile(_src(c), f"{nb.name}::cell{i}", "exec"), ns)
    finally:
        os.chdir(cwd)
    np.save(cache, np.asarray(ns["h_gw_d"], dtype=float))
    return cache


# ===========================================================================
# Subprocess trial runner
# ===========================================================================


def _run_trial(
    soil_tag: str, setup: str, solver: str, model_params: dict, solver_params: dict, et_stress: dict
) -> tuple[dict[str, float], str]:
    """Run a single OAT trial in a subprocess and return all phase metrics.

    Returns a tuple (metrics_dict, error_str).
    metrics_dict has keys: "overall", "wet", "dry", "dry_cool".
    On error, all values are nan and error_str is non-empty.
    """
    pfile = OUT_DIR / "tmp" / f"params-{os.getpid()}-{int(time.time()*1e6)}.json"
    pfile.parent.mkdir(parents=True, exist_ok=True)
    pfile.write_text(
        json.dumps(
            dict(
                soil_tag=soil_tag,
                setup=setup,
                solver=solver,
                model_params=model_params,
                solver_params=solver_params,
                et_stress=et_stress,
            )
        )
    )
    try:
        env = dict(os.environ)
        env.setdefault("MPLCONFIGDIR", str(Path.home() / ".vscode" / "tmp" / "mpl"))
        env.setdefault("KMP_AFFINITY", "disabled")
        proc = subprocess.run(
            [sys.executable, str(WORKER), str(pfile)],
            capture_output=True,
            text=True,
            timeout=WORKER_TIMEOUT_S,
            env=env,
        )
        out = proc.stdout.strip().splitlines()
        for line in out:
            if line.startswith("METRICS_JSON"):
                try:
                    metrics = json.loads(line[len("METRICS_JSON ") :])
                    # Ensure all expected keys present; fill missing with nan.
                    for k in ("overall", "wet", "dry", "dry_cool"):
                        metrics.setdefault(k, float("nan"))
                    return metrics, ""
                except Exception as exc:
                    return dict(_NAN_METRICS), f"METRICS_JSON parse error: {exc}"
            if line.startswith("RMSE_CM"):  # backward compat fallback
                v = float(line.split()[1])
                return {"overall": v, "wet": float("nan"), "dry": float("nan"), "dry_cool": float("nan")}, ""
            if line.startswith("ERROR"):
                return dict(_NAN_METRICS), line
        return dict(_NAN_METRICS), f"no-result rc={proc.returncode} stderr={proc.stderr[-200:]}"
    except subprocess.TimeoutExpired:
        return dict(_NAN_METRICS), "TIMEOUT"
    except Exception as e:
        return dict(_NAN_METRICS), f"{type(e).__name__}: {e}"
    finally:
        try:
            pfile.unlink()
        except Exception:
            pass


# ===========================================================================
# Coordinate descent
# ===========================================================================


def _grids_for(setup: str, solver: str) -> dict:
    g = dict(MODEL_GRID)
    g.update(EXPLICIT_GRID if solver == "explicit" else IMPLICIT_GRID)
    g.update(ET_STRESS_GRID)
    return g


def _accept_tol(setup: str) -> float:
    return IMPROVE_TOL_INT if setup == "intensive" else IMPROVE_TOL_BASIC


def oat_for_case(soil_tag: str, setup: str, solver: str) -> dict:
    print(f"\n=== {soil_tag} / {setup} / {solver} ===", flush=True)
    t_start = time.time()
    _ensure_hydrus_cache(soil_tag, setup)
    base = _parse_baselines(soil_tag, setup)
    mp = dict(base["model_params"])
    et = dict(base["et_stress"])
    sp = dict(base[f"{solver}_solver_params"])
    print(f"  baseline parsed: MP={mp}  ET={et}  SP={sp}", flush=True)

    ref_metrics, err_ref = _run_trial(soil_tag, setup, solver, mp, sp, et)
    rmse_ref = ref_metrics.get("overall", float("nan"))
    if not np.isfinite(rmse_ref):
        return dict(
            soil=soil_tag,
            setup=setup,
            solver=solver,
            error=f"baseline failed: {err_ref}",
            baseline_solver_params=sp,
            baseline_model_params=mp,
            baseline_et_stress=et,
        )
    print(f"  baseline metrics = {ref_metrics}", flush=True)

    grid = _grids_for(setup, solver)
    accept_tol = _accept_tol(setup)

    sweeps: dict = {}
    accepted: list = []
    pass_log: list = []

    # cur_metrics tracks all phase metrics for the current running configuration.
    # Acceptance uses the phase-targeted metric; convergence uses overall.
    cur_metrics = dict(ref_metrics)
    cur_mp, cur_et, cur_sp = dict(mp), dict(et), dict(sp)

    for pass_idx in range(1, MAX_PASSES + 1):
        improved_in_pass = False
        pass_picks = []
        pass_start_overall = cur_metrics["overall"]

        for pname, levels in grid.items():
            in_solver = pname in SOLVER_KEYS
            in_model = pname in cur_mp
            in_et = pname in cur_et
            if not (in_solver or in_model or in_et):
                if pname in MODEL_GRID and pname not in cur_mp:
                    cur_mp[pname] = MODEL_GRID[pname][len(MODEL_GRID[pname]) // 2]
                    in_model = True
                else:
                    continue

            # Determine which metric this parameter should be scored against.
            metric_key = _param_metric(solver, pname)
            cur_metric_val = cur_metrics.get(metric_key, cur_metrics["overall"])
            if not np.isfinite(cur_metric_val):
                # Phase metric unavailable; fall back to overall.
                metric_key = "overall"
                cur_metric_val = cur_metrics["overall"]

            sweep = []
            for lvl in levels:
                mp2, et2, sp2 = dict(cur_mp), dict(cur_et), dict(cur_sp)
                if in_solver:
                    sp2[pname] = lvl
                    if pname == "picard_tol" and lvl < 1e-5:
                        sp2["picard_max_iter"] = 300
                elif in_model:
                    mp2[pname] = lvl
                elif in_et:
                    et2[pname] = lvl
                t0 = time.time()
                trial_metrics, err = _run_trial(soil_tag, setup, solver, mp2, sp2, et2)
                r_targeted = trial_metrics.get(metric_key, float("nan"))
                r_overall = trial_metrics.get("overall", float("nan"))
                ok = np.isfinite(r_targeted)
                sweep.append(
                    dict(
                        value=lvl,
                        rmse=r_targeted,  # phase-targeted metric for this parameter
                        rmse_overall=r_overall,
                        metrics=trial_metrics,
                        metric_key=metric_key,
                        ok=ok,
                        err=err,
                        t_s=time.time() - t0,
                    )
                )

            ok_sweep = [s for s in sweep if s["ok"]]
            if not ok_sweep:
                pass_picks.append((pname, None, None, 0.0))
                continue
            best = min(ok_sweep, key=lambda s: s["rmse"])
            improvement = (cur_metric_val - best["rmse"]) / cur_metric_val if cur_metric_val > 0 else 0.0
            sweeps[pname] = dict(
                levels=sweep,
                best=best,
                pass_idx=pass_idx,
                improvement=improvement,
                metric_key=metric_key,
            )
            if best["rmse"] < cur_metric_val * (1.0 - accept_tol):
                if in_solver:
                    cur_sp[pname] = best["value"]
                    if pname == "picard_tol" and best["value"] < 1e-5:
                        cur_sp["picard_max_iter"] = 300
                elif in_model:
                    cur_mp[pname] = best["value"]
                elif in_et:
                    cur_et[pname] = best["value"]
                # Update all phase metrics for the new running configuration.
                cur_metrics = dict(best["metrics"])
                accepted.append((pname, best["value"], best["rmse"], metric_key))
                pass_picks.append((pname, best["value"], best["rmse"], improvement, metric_key))
                print(
                    f"  [P{pass_idx}] {pname}={best['value']!r}: "
                    f"{metric_key}={best['rmse']:.3f} "
                    f"overall={best['rmse_overall']:.3f} "
                    f"(\u0394={improvement:.1%}) ACCEPTED",
                    flush=True,
                )
            else:
                pass_picks.append((pname, best["value"], best["rmse"], improvement, metric_key))

        # Convergence: check if overall RMSE improved meaningfully this pass.
        if (pass_start_overall - cur_metrics["overall"]) / pass_start_overall >= CONVERGE_TOL:
            improved_in_pass = True

        pass_log.append(pass_picks)
        print(f"  pass {pass_idx} complete -- current overall RMSE = {cur_metrics['overall']:.3f} cm", flush=True)
        if not improved_in_pass:
            print(f"  converged after {pass_idx} pass(es).", flush=True)
            break

    return dict(
        soil=soil_tag,
        setup=setup,
        solver=solver,
        baseline_rmse=rmse_ref,
        baseline_metrics=ref_metrics,
        baseline_solver_params=sp,
        baseline_model_params=mp,
        baseline_et_stress=et,
        final_rmse=cur_metrics["overall"],
        final_metrics=cur_metrics,
        final_solver_params=cur_sp,
        final_model_params=cur_mp,
        final_et_stress=cur_et,
        accepted=accepted,
        pass_log=pass_log,
        sweeps=sweeps,
        wall_s=time.time() - t_start,
    )


# ===========================================================================
# Driver
# ===========================================================================


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="")
    ap.add_argument("--setups", default="basic,intensive")
    ap.add_argument("--solvers", default="implicit,explicit")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--out", default=str(OUT_DIR / "oat_results.json"))
    args = ap.parse_args()

    cases = ALL_CASES if (args.all or not args.cases) else args.cases.split(",")
    setups = args.setups.split(",")
    solvers = args.solvers.split(",")

    results = []
    for soil in cases:
        for setup in setups:
            for solver in solvers:
                try:
                    r = oat_for_case(soil, setup, solver)
                except Exception as e:
                    print(f"FATAL ({soil}/{setup}/{solver}): {e}", flush=True)
                    r = dict(soil=soil, setup=setup, solver=solver, fatal=str(e))
                results.append(r)
                Path(args.out).write_text(json.dumps(results, indent=2, default=str))

    print(f"\nWrote {len(results)} OAT results to {args.out}")


if __name__ == "__main__":
    main()

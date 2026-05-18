"""Apply OAT (iterated coordinate-descent) optima to gen_variants config.

Reads:  tests/sensitivity-analysis/oat_results/oat_results.json
Writes: tests/sensitivity-analysis/oat_results/tuned_params.json

Schema produced by the new oat_harness.py:
    {
      "soil": ..., "setup": ..., "solver": ...,
      "baseline_rmse": ...,   "final_rmse": ...,
      "baseline_metrics": {}, "final_metrics": {},   <- per-phase RMSE dicts
      "baseline_solver_params": {...},  "final_solver_params": {...},
      "baseline_model_params":  {...},  "final_model_params":  {...},
      "baseline_et_stress":     {...},  "final_et_stress":     {...},
      "accepted": [(pname, value, rmse, metric_key), ...],
      "pass_log": [...], "sweeps": {...},
    }

The oat_harness already enforces per-setup acceptance margins
(IMPROVE_TOL_BASIC = 2 %, IMPROVE_TOL_INT = 5 %) during the descent, so
this script simply takes the harness's `final_*` dicts as the chosen
config and confirms there is at least a meaningful net improvement
relative to baseline (using overall RMSE as the gating criterion).
Per-phase metrics are carried through into tuned_params.json for
downstream reporting.
"""

from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "oat_results.json"
OUT = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "tuned_params.json"

NET_IMPROVE_BASIC = 0.02
NET_IMPROVE_INT = 0.05


def pick(case: dict) -> dict:
    setup = case["setup"]
    base = case.get("baseline_rmse")
    final = case.get("final_rmse")
    threshold = NET_IMPROVE_INT if setup == "intensive" else NET_IMPROVE_BASIC

    sd_b = dict(case.get("baseline_solver_params", {}))
    md_b = dict(case.get("baseline_model_params", {}))
    et_b = dict(case.get("baseline_et_stress", {}))
    sd = dict(case.get("final_solver_params", sd_b))
    md = dict(case.get("final_model_params", md_b))
    et = dict(case.get("final_et_stress", et_b))

    note = ""
    if base is None or final is None:
        note = "missing baseline/final; defaults retained"
        sd, md, et = sd_b, md_b, et_b
    else:
        net_improve = (base - final) / base if base > 0 else 0.0
        if net_improve < threshold:
            note = f"net improvement {net_improve:.2%} < threshold {threshold:.0%}; " f"reverting to baseline"
            sd, md, et = sd_b, md_b, et_b
        else:
            note = f"accepted: net {net_improve:.2%} >= {threshold:.0%}"

    return dict(
        solver_params=sd,
        model_params=md,
        et_stress=et,
        baseline_rmse=base,
        final_rmse=final,
        baseline_metrics=case.get("baseline_metrics", {}),
        final_metrics=case.get("final_metrics", {}),
        accepted=case.get("accepted", []),
        note=note,
    )


def main() -> None:
    if not RES.exists():
        raise SystemExit(f"OAT results not found: {RES}")
    rows = json.loads(RES.read_text())

    tuned: dict = {}
    for r in rows:
        if "fatal" in r or "error" in r:
            print(f"SKIP {r['soil']}/{r['setup']}/{r['solver']}: {r.get('error') or r.get('fatal')}")
            continue
        soil, setup, solver = r["soil"], r["setup"], r["solver"]
        cfg = pick(r)
        tuned.setdefault(soil, {}).setdefault(setup, {})[solver] = cfg
        b = cfg.get("baseline_rmse")
        f = cfg.get("final_rmse")
        bs = f"{b:.2f}" if b is not None else "-"
        fs = f"{f:.2f}" if f is not None else "-"
        print(f"{soil:11s} {setup:9s} {solver:8s}  base={bs}cm  final={fs}cm  {cfg['note']}")

    OUT.write_text(json.dumps(tuned, indent=2, default=str))
    print(f"\nWrote tuned config to {OUT}")


if __name__ == "__main__":
    main()

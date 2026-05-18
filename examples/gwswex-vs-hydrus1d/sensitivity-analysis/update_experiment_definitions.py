"""Update experiment_definitions.json with OAT-tuned parameter summary.

Reads:  oat_results/tuned_params.json        (output of apply_oat_optima.py)
Writes: examples/gwswex-vs-hydrus1d/experiment_definitions.json
        (adds / replaces the top-level "oat_tuned" key)

The "oat_tuned" section records the per-(soil, setup, solver) calibrated
configuration for human reference and for cross-validation. It does NOT
replace the fixed geometry / schedule / soil-property / forcing entries,
which remain the single source of truth for the experimental design.

The pipeline writes to tuned_params.json are the authoritative machine-
readable form consumed by gen_variants.py; the oat_tuned section here is
a human-readable record of the same information, co-located with the
fixed experiment definitions for convenience.

Schema added to experiment_definitions.json:
  "oat_tuned": {
    "<soil>": {
      "<setup>": {
        "<solver>": {
          "model_params":  {...},
          "et_stress":     {...},
          "solver_params": {...},
          "baseline_rmse": float,
          "final_rmse":    float,
          "baseline_metrics": {"overall": float, "wet": float, "dry": float, "dry_cool": float},
          "final_metrics":    {"overall": float, "wet": float, "dry": float, "dry_cool": float},
          "note": str,
        }
      }
    }
  }
"""

from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TUNED = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "tuned_params.json"
EXDEF = ROOT / "examples" / "gwswex-vs-hydrus1d" / "experiment_definitions.json"


def main() -> None:
    if not TUNED.exists():
        raise SystemExit(f"tuned_params.json not found: {TUNED}\nRun apply_oat_optima.py first.")

    tuned: dict = json.loads(TUNED.read_text())
    exdef: dict = json.loads(EXDEF.read_text())

    oat_section: dict = {}
    for soil, setups in tuned.items():
        oat_section[soil] = {}
        for setup, solvers in setups.items():
            oat_section[soil][setup] = {}
            for solver, cfg in solvers.items():
                oat_section[soil][setup][solver] = dict(
                    model_params=cfg.get("model_params", {}),
                    et_stress=cfg.get("et_stress", {}),
                    solver_params=cfg.get("solver_params", {}),
                    baseline_rmse=cfg.get("baseline_rmse"),
                    final_rmse=cfg.get("final_rmse"),
                    baseline_metrics=cfg.get("baseline_metrics", {}),
                    final_metrics=cfg.get("final_metrics", {}),
                    note=cfg.get("note", ""),
                )
                b = cfg.get("baseline_rmse"); f = cfg.get("final_rmse")
                bs = f"{b:.3f}" if b is not None else "-"
                fs = f"{f:.3f}" if f is not None else "-"
                print(f"  {soil:11s} {setup:9s} {solver:8s}  {bs} → {fs} cm  {cfg.get('note','')}")

    exdef["oat_tuned"] = oat_section
    exdef["_notes"] = list(exdef.get("_notes", [])) + [
        "oat_tuned: per-(soil, setup, solver) calibrated parameter configurations from "
        "phase-targeted OAT coordinate-descent (sensitivity-analysis/oat_harness.py). "
        "Machine-readable form is in oat_results/tuned_params.json; this section is a "
        "human-readable reference. Fixed geometry / schedule / soil / forcing entries "
        "are unchanged."
    ]
    # Deduplicate _notes
    seen = set()
    exdef["_notes"] = [n for n in exdef["_notes"] if not (n in seen or seen.add(n))]

    EXDEF.write_text(json.dumps(exdef, indent=2) + "\n")
    print(f"\nUpdated {EXDEF.name} with oat_tuned section "
          f"({sum(len(so) for s in oat_section.values() for so in s.values())} entries).")


if __name__ == "__main__":
    main()

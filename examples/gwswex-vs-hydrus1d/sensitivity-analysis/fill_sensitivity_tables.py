"""Fill §2.3 per-case parameter tables in
docs/sensitivity-analysis-from-comparison.md from tuned_params.json.

Reads:  tests/sensitivity-analysis/oat_results/tuned_params.json
Writes: docs/sensitivity-analysis-from-comparison.md (in-place)

The script replaces the placeholder rows in §2.3.1 (implicit) and
§2.3.2 (explicit) and leaves the rest of the doc untouched.
"""

from __future__ import annotations
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TUNED = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "tuned_params.json"
DOC = ROOT / "docs" / "sensitivity-analysis-from-comparison.md"

SOILS = ["loam", "sand", "clay", "sand-loam", "sand-clay", "loam-clay"]
SETUPS = ["basic", "intensive"]

IMPL_COLS = ("picard_tol", "picard_max_iter", "n_trapz", "beta_hyst", "psi_f", "F_min", "ICratio_min")
EXPL_COLS = ("courant_number", "n_trapz", "beta_hyst", "psi_f", "F_min", "ICratio_min")
ET_COLS = ("s_star", "s_w", "s_h", "s_e")


def _fmt(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e4):
            return f"{v:.0e}"
        return f"{v:g}"
    return str(v)


def _row(soil: str, setup: str, cfg: dict, cols) -> str:
    sp = cfg.get("solver_params", {}) or {}
    mp = cfg.get("model_params", {}) or {}
    et = cfg.get("et_stress", {}) or {}
    rmse = cfg.get("final_rmse")
    cells = [soil, setup]
    for k in cols:
        if k in sp:
            cells.append(_fmt(sp[k]))
        elif k in mp:
            cells.append(_fmt(mp[k]))
        else:
            cells.append("—")
    for k in ET_COLS:
        cells.append(_fmt(et.get(k)))
    cells.append(f"{rmse:.2f}" if rmse is not None else "—")
    return "| " + " | ".join(cells) + " |"


def main() -> None:
    if not TUNED.exists():
        raise SystemExit(f"tuned_params.json not found at {TUNED}")
    tuned = json.loads(TUNED.read_text())

    impl_rows = []
    expl_rows = []
    for soil in SOILS:
        for setup in SETUPS:
            cell = tuned.get(soil, {}).get(setup, {})
            if "implicit" in cell:
                impl_rows.append(_row(soil, setup, cell["implicit"], IMPL_COLS))
            if "explicit" in cell:
                expl_rows.append(_row(soil, setup, cell["explicit"], EXPL_COLS))

    impl_table = (
        "#### 2.3.1 Implicit solver\n\n"
        "| Soil | Setup | `picard_tol` | `picard_max_iter` | `n_trapz` | `beta_hyst` | "
        "`psi_f` | `F_min` | `ICratio_min` | $s^*$ | $s_w$ | $s_h$ | $s_e$ | RMSE [cm] |\n"
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n" + "\n".join(impl_rows) + "\n"
    )
    expl_table = (
        "#### 2.3.2 Explicit solver\n\n"
        "| Soil | Setup | `courant_number` | `n_trapz` | `beta_hyst` | "
        "`psi_f` | `F_min` | `ICratio_min` | $s^*$ | $s_w$ | $s_h$ | $s_e$ | RMSE [cm] |\n"
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|\n" + "\n".join(expl_rows) + "\n"
    )

    doc = DOC.read_text()
    doc = re.sub(
        r"#### 2\.3\.1 Implicit solver[\s\S]+?(?=#### 2\.3\.2)",
        impl_table + "\n",
        doc,
    )
    doc = re.sub(
        r"#### 2\.3\.2 Explicit solver[\s\S]+?(?=\n## 3\. Discussion)",
        expl_table + "\n",
        doc,
    )
    DOC.write_text(doc)
    print(f"Wrote {len(impl_rows)} implicit + {len(expl_rows)} explicit rows to {DOC.name}")


if __name__ == "__main__":
    main()

"""Validate all 12 comparison notebooks: run end-to-end, scrape RMSE/NSE.

Writes summary to tests/sensitivity-analysis/oat_results/validation_summary.md

Metric lines produced by the executed notebooks:
  WT depth     Full (0-65 d)    implicit       0.63 cm     0.42 cm    0.999    +0.42 cm
  WT depth     Wet (5-35 d)     implicit       0.73 cm     0.51 cm    0.996    +0.51 cm
  WT depth     Dry (35-65 d)    implicit       0.57 cm     0.40 cm    0.962    +0.40 cm
"""

from __future__ import annotations
import json
import re
import subprocess
import time
from pathlib import Path

import os
import sys

ROOT = Path(__file__).resolve().parents[3]
EX = ROOT / "examples" / "gwswex-vs-hydrus1d"
LOG_DIR = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "validation_logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "validation_summary.md"

SOILS = ["loam", "sand", "clay", "sand-loam", "sand-clay", "loam-clay"]
SETUPS = ["basic", "intensive"]

# Regex patterns for Full, Wet, and Dry metric lines.
_METRIC_RE = re.compile(
    r"WT depth\s+(Full|Wet|Dry)\s*\([^)]*\)\s+(implicit|explicit)\s+" r"([\d.]+)\s*cm\s+[\d.]+\s*cm\s+(-?[\d.]+)"
)
# Backward-compat alias used below for full-horizon scraping.
LINE_RE = _METRIC_RE

rows = []
for setup in SETUPS:
    for soil in SOILS:
        nb = EX / (f"comparison-basic-{soil}.ipynb" if setup == "basic" else f"comparison-intensive-{soil}.ipynb")
        if not nb.exists() and setup == "basic":
            nb = EX / f"comparison-{soil}.ipynb"
        log = LOG_DIR / f"val-{setup}-{soil}.log"
        t0 = time.time()
        rc = subprocess.call(
            [
                sys.executable,
                "-m",
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--inplace",
                "--ExecutePreprocessor.timeout=600",
                str(nb),
            ],
            stdout=log.open("w"),
            stderr=subprocess.STDOUT,
            cwd=str(EX),
        )
        dt = time.time() - t0
        # Parse executed notebook JSON for stdout outputs (nbconvert --inplace
        # writes them into the notebook, not the subprocess stdout).
        try:
            nb_json = json.loads(nb.read_text())
            stream_text = []
            for cell in nb_json.get("cells", []):
                for out in cell.get("outputs", []) or []:
                    if out.get("output_type") == "stream":
                        s = out.get("text", "")
                        stream_text.append("".join(s) if isinstance(s, list) else s)
            text = "\n".join(stream_text)
        except Exception:
            text = log.read_text(errors="ignore")
        # Collect per-phase blocks: keyed by (phase_label, solver).
        phase_data: dict[tuple[str, str], tuple[float, float]] = {}
        for ln in text.splitlines():
            m = _METRIC_RE.search(ln)
            if m:
                phase_label = m.group(1).lower()  # "full", "wet", "dry"
                solver = m.group(2)  # "implicit", "explicit"
                rmse_val = float(m.group(3))
                nse_val = float(m.group(4))
                phase_data[(phase_label, solver)] = (rmse_val, nse_val)

        def _get(phase: str, solver: str) -> tuple:
            return phase_data.get((phase, solver), (None, None))

        rows.append(
            dict(
                setup=setup,
                soil=soil,
                rc=rc,
                dt=dt,
                impl_rmse=_get("full", "implicit")[0],
                impl_nse=_get("full", "implicit")[1],
                expl_rmse=_get("full", "explicit")[0],
                expl_nse=_get("full", "explicit")[1],
                impl_wet_rmse=_get("wet", "implicit")[0],
                impl_wet_nse=_get("wet", "implicit")[1],
                expl_wet_rmse=_get("wet", "explicit")[0],
                expl_wet_nse=_get("wet", "explicit")[1],
                impl_dry_rmse=_get("dry", "implicit")[0],
                impl_dry_nse=_get("dry", "implicit")[1],
                expl_dry_rmse=_get("dry", "explicit")[0],
                expl_dry_nse=_get("dry", "explicit")[1],
            )
        )
        print(
            f"{setup:9s} {soil:11s} rc={rc} {dt:6.1f}s   "
            f"impl={_get('full','implicit')[0] or 'NA'}cm/{_get('full','implicit')[1] or 'NA'}   "
            f"expl={_get('full','explicit')[0] or 'NA'}cm/{_get('full','explicit')[1] or 'NA'}",
            flush=True,
        )


def _fmt(v) -> str:
    return f"{v:.2f}" if v is not None else "FAIL"


# Markdown summary
with OUT.open("w") as f:
    f.write("# Validation summary — OAT-tuned comparison notebooks\n\n")
    f.write("## Full-horizon RMSE and NSE\n\n")
    f.write("| setup | soil | rc | runtime [s] | impl RMSE [cm] | impl NSE | expl RMSE [cm] | expl NSE |\n")
    f.write("|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        f.write(
            f"| {r['setup']} | {r['soil']} | {r['rc']} | {r['dt']:.1f} | "
            f"{_fmt(r['impl_rmse'])} | {_fmt(r['impl_nse'])} | "
            f"{_fmt(r['expl_rmse'])} | {_fmt(r['expl_nse'])} |\n"
        )
    f.write("\n## Phase-resolved RMSE and NSE\n\n")
    f.write("| setup | soil | solver | wet RMSE [cm] | wet NSE | dry RMSE [cm] | dry NSE |\n")
    f.write("|---|---|---|---|---|---|---|\n")
    for r in rows:
        for slv in ("impl", "expl"):
            f.write(
                f"| {r['setup']} | {r['soil']} | {slv} | "
                f"{_fmt(r[f'{slv}_wet_rmse'])} | {_fmt(r[f'{slv}_wet_nse'])} | "
                f"{_fmt(r[f'{slv}_dry_rmse'])} | {_fmt(r[f'{slv}_dry_nse'])} |\n"
            )
print(f"\nWrote {OUT}")

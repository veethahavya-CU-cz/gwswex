"""Single-trial GWSWEX worker for OAT harness.

Runs in a fresh subprocess (one per trial) to isolate the global Fortran
kernel state. If a trial crashes, only this subprocess dies; the parent
harness reads `nan` and moves on.

Reads JSON from argv[1] (path to params file). Writes two lines to stdout:
  METRICS_JSON {"overall": x, "wet": y, "dry": z, "dry_cool": w}
  RMSE_CM <overall>  (backward-compat line for harness fallback)
Or on error:
  ERROR <msg>

Phase boundaries (0-based indices into the output array):
  basic (65-element daily array):
    wet  = slice(5, 35)   → days 6–35  (T_P2=5 d, T_P3=35 d)
    dry  = slice(35, 65)  → days 36–65
    dry_cool = dry        (no cooldown phase in basic)
  intensive (768-element hourly array):
    wet  = slice(72, 312)  → hours 72–311  (T_WU=72 h, T_WET=240 h)
    dry  = slice(312, 600) → hours 312–599 (T_DRY=288 h)
    cool = slice(600, 768) → hours 600–767 (T_COOL=168 h)
    dry_cool = slice(312, 768)

NaN-aware RMSE: HYDRUS-1D returns NaN during ponding (intensive setup).
All phase metrics use _rmse_masked() which excludes NaN-contaminated steps.

Params JSON schema:
    {
      "soil_tag": "loam",
      "setup":    "basic" | "intensive",
      "solver":   "implicit" | "explicit",
      "model_params": {...},
      "solver_params": {...},
      "et_stress":   {...},
    }
"""

from __future__ import annotations

import json
import os
import sys
import traceback
from datetime import timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.show = lambda: None  # type: ignore

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
EX_DIR = ROOT / "examples" / "gwswex-vs-hydrus1d"


def _nb_name(soil_tag: str, setup: str) -> str:
    if setup == "basic":
        nb = EX_DIR / f"comparison-basic-{soil_tag}.ipynb"
        if not nb.exists():
            nb = EX_DIR / f"comparison-{soil_tag}.ipynb"
        return nb.name
    return f"comparison-intensive-{soil_tag}.ipynb"


CACHE_DIR = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "cache"


def _build_namespace(soil_tag: str, setup: str) -> dict:
    """Exec cells 0, 1, 3 of the notebook; replace cell 2 (HYDRUS) with a
    cached `h_gw_d` load so the worker is fast (~1 s vs ~10 s)."""
    nb = EX_DIR / _nb_name(soil_tag, setup)
    cells = json.loads(nb.read_text())["cells"]
    code_cells = [c for c in cells if c.get("cell_type") == "code"]
    ns: dict = {"__name__": "__oat_worker__", "__file__": str(nb), "__setup__": setup}
    cache = CACHE_DIR / f"h_gw_d-{soil_tag}-{setup}.npy"
    cwd = os.getcwd()
    os.chdir(EX_DIR)
    try:
        for i, c in enumerate(code_cells[:4]):
            src = c["source"]
            if isinstance(src, list):
                src = "".join(src)
            if i == 2 and cache.exists():
                src = "import numpy as _np_oat\n" f"h_gw_d = _np_oat.load(r'{cache}')\n"
            exec(compile(src, f"{nb.name}::cell{i}", "exec"), ns)
    finally:
        os.chdir(cwd)
    return ns


def _hydrus_ref(ns: dict) -> np.ndarray:
    # Defensive clip: ensure HYDRUS WT depth is pinned at the surface (>=0 cm)
    # before scoring. The notebooks do this in cell 2 right after np.interp,
    # but we re-apply here so the harness scoring is unambiguous and matches
    # the notebook ground truth even if the cache was built before clipping
    # was added or if a future notebook edit drops it.
    return np.maximum(np.asarray(ns["h_gw_d"], dtype=float), 0.0)


def _rmse_masked(a: np.ndarray, b: np.ndarray) -> float:
    """NaN-aware RMSE. Excludes timesteps where either a or b is non-finite.
    Returns nan when no valid timesteps remain (e.g. empty phase slice)."""
    n = min(len(a), len(b))
    diff = (a[:n] - b[:n]).astype(float)
    valid = np.isfinite(diff)
    if not valid.any():
        return float("nan")
    return float(np.sqrt(np.mean(diff[valid] ** 2)))


# Keep legacy _rmse for internal use only; all public scoring uses _rmse_masked.
def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    n = min(len(a), len(b))
    return float(np.sqrt(np.mean((a[:n] - b[:n]) ** 2)))


def _phase_bounds(setup: str) -> dict:
    """Return 0-based slice objects for wet / dry / cool / dry_cool phases.

    Indices into the GWSWEX output array (basic: 65-element daily;
    intensive: 768-element hourly). Derived from experiment_definitions.json:
      basic:     T_P2=5 d (warmup end / wet start), T_P3=35 d (wet end / dry start)
      intensive: T_WU=72 h, T_WET=240 h, T_DRY=288 h, T_COOL=168 h
    """
    if setup == "basic":
        wet = slice(5, 35)  # days 6–35  (0-indexed: 5..34)
        dry = slice(35, 65)  # days 36–65 (0-indexed: 35..64)
        return {"wet": wet, "dry": dry, "cool": None, "dry_cool": dry}
    # intensive
    wet = slice(72, 312)  # hours 72–311
    dry = slice(312, 600)  # hours 312–599
    cool = slice(600, 768)  # hours 600–767
    dry_cool = slice(312, 768)  # dry + cooldown combined
    return {"wet": wet, "dry": dry, "cool": cool, "dry_cool": dry_cool}


def _run(soil_tag: str, setup: str, solver: str, model_params: dict, solver_params: dict, et_stress: dict) -> dict:
    from gwswex.model import GWSWEXmodel  # type: ignore

    ns = _build_namespace(soil_tag, setup)
    Z_TOP = ns["Z_TOP"]
    Z_WT = ns["Z_WT"]
    Z_BOT = ns["Z_BOT"]
    DZ = ns["DZ"]
    NL = ns["NL"]
    LAYER_MAT_IDS = ns["LAYER_MAT_IDS"]
    MATERIALS = ns["MATERIALS"]
    ROOT_D0 = ns["ROOT_D0"]
    ROOT_D1 = ns["ROOT_D1"]

    if setup == "intensive":
        N = ns["N"]
        T_unit = "h"
        bot = [list(np.round(np.linspace(Z_TOP - DZ, Z_BOT, NL), 10))]
        mats_local = []
        for _mat in MATERIALS:
            _m = dict(_mat)
            _m["K_sat"] = _mat["K_sat"] / 24.0
            mats_local.append(_m)
        forcing = None
    else:
        BNDS = ns["BNDS"]
        T_TOTAL = ns["T_TOTAL"]
        T0 = ns["T0"]
        forcing = ns["f_e"] if solver == "explicit" else ns["f_i"]
        T_unit = "d"
        bot = [list(BNDS[1:])]
        mats_local = MATERIALS

    out = (
        ROOT
        / "examples"
        / "gwswex-vs-hydrus1d"
        / "sensitivity-analysis"
        / "oat_results"
        / "tmp"
        / f"oat-{os.getpid()}.nc"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()

    m = GWSWEXmodel(name=f"oat-{solver}", T=T_unit, L="m", output_fpath=str(out))
    m.init_space(ne=1, nl=NL, top=[[Z_TOP]], bot=bot, sID=[LAYER_MAT_IDS], vID=[[1]])
    for _mat in mats_local:
        m.add_material(**_mat)
    m.add_vegetation(
        id=1,
        name="crop",
        root_depth_initial=ROOT_D0,
        root_depth_final=ROOT_D1,
        root_growth_model="linear",
        et_stress=et_stress,
    )
    if setup == "intensive":
        m.init_time(n_steps=N, dt=1.0, dt_min=1 / 60)
    else:
        m.init_time(
            start=T0,
            stop=T0 + timedelta(days=T_TOTAL),
            dt=timedelta(hours=1) if solver == "implicit" else timedelta(days=1),
            dt_min=timedelta(seconds=60),
            adaptive=True,
        )
    m.set_model_params(**model_params)
    m.set_solver(solver=solver, **solver_params)
    m.set_initial_conditions(
        gw=[Z_WT] if setup == "intensive" else Z_WT,
        sw=[0.0] if setup == "intensive" else 0.0,
        uz=[-999] if setup == "intensive" else -999,
    )
    if isinstance(forcing, dict):
        m.set_forcing(**forcing)
    else:
        m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
    m.init()

    gw_all: list = []
    h_ref = _hydrus_ref(ns)
    if setup == "intensive":
        prec_h = ns["prec_h"]
        pet_h = ns["pet_h"]
        ptt_h = ns["ptt_h"]
        for t in range(N):
            m.update_forcing(t, precip=float(prec_h[t]), pet=float(pet_h[t]), ptt=float(ptt_h[t]))
            m.run_step(t)
            gw_all.append(float(m.get_state()["GWH"][0]))
        m.deinit()
        gw_arr = np.array(gw_all)
        wt_clipped = np.maximum((Z_TOP - gw_arr) * 100.0, 0.0)
    else:

        def _cb(t, state):
            gw_all.append(state["GWH"][0])

        m.run(callback=_cb)
        m.deinit()
        gw_arr = np.array(gw_all)
        steps_per_day = 24 if solver == "implicit" else 1
        idx = np.arange(steps_per_day - 1, len(gw_arr), steps_per_day)
        wt_clipped = np.maximum((Z_TOP - gw_arr[idx]) * 100.0, 0.0)

    bounds = _phase_bounds(setup)
    overall = _rmse_masked(wt_clipped, h_ref)
    wet = _rmse_masked(wt_clipped[bounds["wet"]], h_ref[bounds["wet"]])
    dry = _rmse_masked(wt_clipped[bounds["dry"]], h_ref[bounds["dry"]])
    dry_cool = _rmse_masked(wt_clipped[bounds["dry_cool"]], h_ref[bounds["dry_cool"]])
    return {"overall": overall, "wet": wet, "dry": dry, "dry_cool": dry_cool}


def main() -> None:
    if len(sys.argv) != 2:
        print("ERROR usage: oat_worker.py <params.json>")
        sys.exit(2)
    p = json.loads(Path(sys.argv[1]).read_text())
    try:
        metrics = _run(
            soil_tag=p["soil_tag"],
            setup=p["setup"],
            solver=p["solver"],
            model_params=p["model_params"],
            solver_params=p["solver_params"],
            et_stress=p["et_stress"],
        )
        print(f"METRICS_JSON {json.dumps(metrics)}")
        print(f"RMSE_CM {metrics['overall']:.6f}")  # backward-compat
    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        print(f"ERROR {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()

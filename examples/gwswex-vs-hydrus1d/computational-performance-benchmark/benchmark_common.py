# type: ignore
"""Shared setup for the computational-performance benchmark.

Loads the single source of truth (`experiment_definitions.json`) and derives
the intensive-sand-loam configuration used by both `run_gwswex.py` and
`run_hydrus.py`. Keeping the loader in one module guarantees that all three
codes (GWSWEX explicit, GWSWEX implicit, HYDRUS-1D) see identical soil,
geometry, schedule and forcing.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EX_DIR = HERE.parent  # examples/gwswex-vs-hydrus1d/
ROOT = EX_DIR.parent.parent  # repository root
RESULTS_DIR = HERE / "results"
TMP_DIR = HERE / "tmp"
SOT_PATH = EX_DIR / "experiment_definitions.json"
HYDRUS_EXE = EX_DIR / "hydrus1d" / "bin" / "hydrus"

SETUP = "intensive"
# SOIL_TAG selects a variant from experiment_definitions.json["variants"].
# Override at runtime via env var GWSWEX_BENCH_SOIL (e.g. "loam", "loam-clay").
SOIL_TAG = os.environ.get("GWSWEX_BENCH_SOIL", "sand-loam")
NE_DEFAULT = (1, 10, 100, 1000, 10000)
N_TRIALS_DEFAULT = 3


def load_sot() -> dict:
    return json.loads(SOT_PATH.read_text())


def derive_intensive_case() -> dict:
    """Return a dict carrying every quantity the benchmark drivers need.

    Time unit: hours. Length unit: metres. The intensive setup is hourly.
    Soil and variant are selected by the module-level `SOIL_TAG`.
    """
    sot = load_sot()
    g = sot["geometry"][SETUP]
    s = sot["schedule"][SETUP]
    veg = sot["vegetation"][SETUP]
    var = sot["variants"][SOIL_TAG]
    layers_spec = var["layers"]  # list of [soil_name, fraction]
    # Primary soil = first listed (used as the reference for `cfg['soil']`).
    primary_name = layers_spec[0][0]
    soil = {k: v for k, v in sot["soils"][primary_name].items() if not k.startswith("_")}

    Z_TOP = float(g["Z_TOP"])
    Z_BOT = float(g["Z_BOT"])
    DZ = float(g["DZ"])
    NL = int(g["NL"])
    Z_WT = float(g["Z_WT"])

    T_WU = int(s["T_WU_h"])
    T_WET = int(s["T_WET_h"])
    T_DRY = int(s["T_DRY_h"])
    T_COOL = int(s["T_COOL_h"])
    N = T_WU + T_WET + T_DRY + T_COOL  # total hourly steps (768 h = 32 d)
    T_TOTAL_D = N // 24

    f = var[SETUP]
    h_idx = np.arange(N)
    prec = np.where(h_idx < T_WU, 0.0, np.where(h_idx < T_WU + T_WET, f["P_W"], 0.0)).astype(float)
    pet = np.where(
        h_idx < T_WU,
        0.0,
        np.where(h_idx < T_WU + T_WET, f["E_W"], np.where(h_idx < T_WU + T_WET + T_DRY, f["E_D"], 0.0)),
    ).astype(float)
    ptt = np.where(
        h_idx < T_WU,
        0.0,
        np.where(h_idx < T_WU + T_WET, f["T_W"], np.where(h_idx < T_WU + T_WET + T_DRY, f["T_D"], 0.0)),
    ).astype(float)

    # Build materials list (one per distinct soil in layers_spec, preserving order)
    materials: list[dict] = []
    name_to_id: dict[str, int] = {}
    for soil_name, _frac in layers_spec:
        if soil_name in name_to_id:
            continue
        s_props = sot["soils"][soil_name]
        mid = len(materials) + 1
        name_to_id[soil_name] = mid
        materials.append(
            dict(
                id=mid,
                name=soil_name,
                K_sat=s_props["K_sat"] / 24.0,  # m d-1 -> m h-1 (intensive is hourly)
                lam=s_props["lam"],
                vanG=dict(
                    alpha=s_props["alpha"], n=s_props["n"], theta_r=s_props["theta_r"], theta_s=s_props["theta_s"]
                ),
            )
        )

    # Layer geometry: assign mat_id per layer based on depth fractions.
    # Layers are top-down (layer 0 is shallowest); fractions are top-down too.
    bot_layer_elevations = list(np.round(np.linspace(Z_TOP - DZ, Z_BOT, NL), 10))
    total_depth = Z_TOP - Z_BOT
    cum_frac = 0.0
    layer_mat_ids: list[int] = []
    boundaries: list[float] = []  # cumulative depth (m) at end of each spec layer
    for _name, frac in layers_spec:
        cum_frac += frac
        boundaries.append(Z_TOP - cum_frac * total_depth)
    # For each numerical layer, assign mat_id by midpoint elevation
    layer_top = Z_TOP
    for i in range(NL):
        layer_bot = bot_layer_elevations[i]
        mid_elev = 0.5 * (layer_top + layer_bot)
        mat_id = name_to_id[layers_spec[-1][0]]  # default = deepest spec layer
        for j, b in enumerate(boundaries):
            if mid_elev >= b - 1e-12:
                mat_id = name_to_id[layers_spec[j][0]]
                break
        layer_mat_ids.append(mat_id)
        layer_top = layer_bot

    # Daily forcing for HYDRUS atmospheric BC (cm d-1)
    M_PER_H_TO_CM_PER_D = 100.0 * 24.0
    day_idx = np.arange(1, T_TOTAL_D + 1, dtype=float)
    P_W_CMD = f["P_W"] * M_PER_H_TO_CM_PER_D
    E_W_CMD = f["E_W"] * M_PER_H_TO_CM_PER_D
    T_W_CMD = f["T_W"] * M_PER_H_TO_CM_PER_D
    E_D_CMD = f["E_D"] * M_PER_H_TO_CM_PER_D
    T_D_CMD = f["T_D"] * M_PER_H_TO_CM_PER_D
    prec_d = np.where((day_idx > 3) & (day_idx <= 13), P_W_CMD, 0.0)
    pet_d = np.where((day_idx > 3) & (day_idx <= 13), E_W_CMD, np.where((day_idx > 13) & (day_idx <= 25), E_D_CMD, 0.0))
    ptt_d = np.where((day_idx > 3) & (day_idx <= 13), T_W_CMD, np.where((day_idx > 13) & (day_idx <= 25), T_D_CMD, 0.0))

    return dict(
        # geometry
        Z_TOP=Z_TOP,
        Z_BOT=Z_BOT,
        DZ=DZ,
        NL=NL,
        Z_WT=Z_WT,
        # schedule
        N=N,
        T_TOTAL_D=T_TOTAL_D,
        T_WU=T_WU,
        T_WET=T_WET,
        T_DRY=T_DRY,
        T_COOL=T_COOL,
        # forcing (hourly, m h-1)
        prec_h=prec,
        pet_h=pet,
        ptt_h=ptt,
        # forcing (daily, cm d-1, for HYDRUS)
        day_idx=day_idx,
        prec_d=prec_d,
        pet_d=pet_d,
        ptt_d=ptt_d,
        # materials
        soil=soil,
        materials=materials,
        layer_mat_ids=layer_mat_ids,
        bot_layer_elevations=bot_layer_elevations,
        # vegetation
        ROOT_D0=float(veg["root_d_initial"]),
        ROOT_D1=float(veg["root_d_final"]),
    )


def ensure_dirs() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

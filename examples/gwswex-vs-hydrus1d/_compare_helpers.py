# type: ignore
"""Shared pre/post/plot helpers for the GWSWEX vs HYDRUS-1D comparison notebooks.

The notebooks themselves are kept focused on GWSWEX model setup and run logic.
Everything else - SoT loading, HYDRUS-1D drive+parse, plotting, metrics - lives
here. Public surface is small and intentionally vanilla-Python so notebook code
reads top-to-bottom without indirection.

Public surface
--------------
load_experiment(setup, soil_tag) -> ExperimentSpec
    Resolves geometry, schedule, soil layering, vegetation, forcing arrays
    and tuned per-solver parameters from experiment_definitions.json.

run_hydrus(spec) -> HydrusResult
    Runs HYDRUS-1D once via phydrus and returns aligned hourly/daily series.

setup_plot_style(), STYLE, PHASES_basic / PHASES_intensive, shade_phases(...)

plot_water_table(spec, hyd, runs) ; plot_zone_theta(spec, hyd, runs, zones)
plot_theta_snapshots(spec, hyd, runs, snaps)
plot_cum_et(spec, hyd, runs) ; plot_surface_ponding(spec, hyd, runs)
plot_mass_balance(spec, hyd, runs) ; plot_metrics_summary(spec, hyd, runs, ...)

rmse, mae, bias, nse  -- scalar metric helpers
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOT_PATH = HERE / "experiment_definitions.json"


# ─────────────────────────────────────────────────────────────────────────────
# Experiment specification (Single source of truth, JSON-driven)
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class ExperimentSpec:
    """Resolved experiment definition for one (setup, soil_tag) pair."""

    setup: str  # "basic" or "intensive"
    soil_tag: str  # variant key (e.g. "sand-loam")

    # Geometry
    Z_TOP: float
    Z_BOT: float
    DZ: float
    NL: int
    Z_WT: float
    BNDS: np.ndarray  # layer interfaces [m], top->bot
    LAYER_DEPTH_CM: np.ndarray  # midpoint depth below surface [cm]

    # Soil layering and materials (top-to-bottom)
    materials: list[dict]  # ready for GWSWEXmodel.add_material(**m)
    soil_db: dict[str, dict]  # all SoT soil records
    layer_mat_ids: list[int]
    theta_r: float  # bottom-material θr (GW zone)
    theta_s: float  # bottom-material θs

    # Schedule
    T_TOTAL_d: int  # total run length [d]
    phase_bounds_d: list[float]  # phase boundary times [d], len=nphases+1
    phase_labels: list[str]
    phase_colors: list[str]

    # Time-axis (hourly grid for intensive, daily for basic)
    t_d: np.ndarray  # plotting axis [d]
    n_steps: int  # GWSWEX macro-steps in run
    dt_unit: str  # GWSWEX time unit ('h' or 'd')

    # Forcing
    prec_native: np.ndarray  # in GWSWEX units (m d^-1 or m h^-1)
    pet_native: np.ndarray
    ptt_native: np.ndarray
    prec_d_cmd: np.ndarray  # daily series for HYDRUS-1D atm BC [cm d^-1]
    pet_d_cmd: np.ndarray
    ptt_d_cmd: np.ndarray
    day_idx: np.ndarray  # 1..T_TOTAL_d

    # Vegetation
    root_d_initial: float
    root_d_final: float
    root_growth_model: str

    # Tuned params per solver (from oat_tuned)
    tuned: dict[str, dict]  # tuned[solver] = {model_params, et_stress, solver_params}

    # Output paths
    outputs_root: Path
    hydrus_ws: Path

    # Reference summary (rendered in markdown headers)
    summary: dict[str, Any] = field(default_factory=dict)


def load_experiment(setup: str, soil_tag: str) -> ExperimentSpec:
    """Build an ExperimentSpec from experiment_definitions.json."""
    sot = json.loads(SOT_PATH.read_text())

    # Geometry --------------------------------------------------------------
    g = sot["geometry"][setup]
    Z_TOP, Z_BOT, DZ = float(g["Z_TOP"]), float(g["Z_BOT"]), float(g["DZ"])
    NL, Z_WT = int(g["NL"]), float(g["Z_WT"])
    BNDS = np.linspace(Z_TOP, Z_BOT, NL + 1)
    LAYER_DEPTH_CM = np.linspace(DZ / 2, (Z_TOP - Z_BOT) - DZ / 2, NL) * 100.0

    # Soil layering ---------------------------------------------------------
    soil_db = {n: {k: v for k, v in d.items() if not k.startswith("_")} for n, d in sot["soils"].items()}
    var = sot["variants"][soil_tag]
    layers = [(name, float(frac)) for name, frac in var["layers"]]

    materials, seen = [], {}
    for name, _ in layers:
        if name in seen:
            continue
        seen[name] = len(materials) + 1
        s = soil_db[name]
        materials.append(
            dict(
                id=seen[name],
                name=name,
                K_sat=s["K_sat"],
                lam=s["lam"],
                vanG=dict(alpha=s["alpha"], n=s["n"], theta_r=s["theta_r"], theta_s=s["theta_s"]),
            )
        )
    layer_mat_ids: list[int] = []
    for name, frac in layers:
        layer_mat_ids += [seen[name]] * int(round(NL * frac))
    layer_mat_ids = (layer_mat_ids + [layer_mat_ids[-1]] * NL)[:NL]
    bot_soil = soil_db[layers[-1][0]]
    theta_r, theta_s = bot_soil["theta_r"], bot_soil["theta_s"]

    # Schedule + forcing ----------------------------------------------------
    sch = sot["schedule"][setup]
    f = var[setup]
    if setup == "basic":
        T_TOTAL_d = int(sch["T_TOTAL_d"])
        T_P1, T_P2, T_P3 = sch["T_P1_d"], sch["T_P2_d"], sch["T_P3_d"]
        phase_bounds_d = [float(T_P1), float(T_P2), float(T_P3), float(T_TOTAL_d)]
        phase_labels = ["Warmup", "Wet", "Dry"]
        phase_colors = ["#ddeef7", "#ceecc8", "#faecd1"]
        t_d = np.arange(1, T_TOTAL_d + 1, dtype=float)
        n_steps, dt_unit = T_TOTAL_d, "d"
        # Daily forcing in cm d^-1 -> GWSWEX m d^-1
        prec_cmd = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, f["P_WET"], 0.0))
        pet_cmd = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, f["PE_WET"], f["PE_DRY"]))
        ptt_cmd = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, f["PT_WET"], f["PT_DRY"]))
        cm2m = 0.01
        prec_native = prec_cmd * cm2m
        pet_native = pet_cmd * cm2m
        ptt_native = ptt_cmd * cm2m
        prec_d_cmd, pet_d_cmd, ptt_d_cmd = prec_cmd, pet_cmd, ptt_cmd
        day_idx = t_d.copy()
        summary = dict(
            P_WET_cmd=f["P_WET"],
            PE_WET_cmd=f["PE_WET"],
            PT_WET_cmd=f["PT_WET"],
            PE_DRY_cmd=f["PE_DRY"],
            PT_DRY_cmd=f["PT_DRY"],
            T_P1=T_P1,
            T_P2=T_P2,
            T_P3=T_P3,
            T_TOTAL=T_TOTAL_d,
        )
    else:  # intensive
        T_WU, T_WET, T_DRY, T_COOL = (
            int(sch["T_WU_h"]),
            int(sch["T_WET_h"]),
            int(sch["T_DRY_h"]),
            int(sch["T_COOL_h"]),
        )
        N = T_WU + T_WET + T_DRY + T_COOL
        T_TOTAL_d = N // 24
        T_P1 = 0.0
        T_P2 = T_WU / 24.0
        T_P3 = (T_WU + T_WET) / 24.0
        T_P4 = (T_WU + T_WET + T_DRY) / 24.0
        phase_bounds_d = [T_P1, T_P2, T_P3, T_P4, float(T_TOTAL_d)]
        phase_labels = ["Warmup", "Wet", "Dry", "Cooldown"]
        phase_colors = ["#ddeef7", "#ceecc8", "#faecd1", "#e8e0f0"]
        t_d = np.arange(1, N + 1) / 24.0
        n_steps, dt_unit = N, "h"
        h_idx = np.arange(N)
        P_W, E_W, T_W = f["P_W"], f["E_W"], f["T_W"]
        E_D, T_D = f["E_D"], f["T_D"]
        prec_native = np.where(h_idx < T_WU, 0.0, np.where(h_idx < T_WU + T_WET, P_W, 0.0))
        pet_native = np.where(
            h_idx < T_WU, 0.0, np.where(h_idx < T_WU + T_WET, E_W, np.where(h_idx < T_WU + T_WET + T_DRY, E_D, 0.0))
        )
        ptt_native = np.where(
            h_idx < T_WU, 0.0, np.where(h_idx < T_WU + T_WET, T_W, np.where(h_idx < T_WU + T_WET + T_DRY, T_D, 0.0))
        )
        # HYDRUS-1D atm BC needs daily totals in cm d^-1
        m_per_h_to_cm_per_d = 100.0 * 24.0
        day_idx = np.arange(1, T_TOTAL_d + 1, dtype=float)
        wet_day_lo, wet_day_hi = T_WU // 24, (T_WU + T_WET) // 24
        dry_day_hi = (T_WU + T_WET + T_DRY) // 24
        prec_d_cmd = np.where((day_idx > wet_day_lo) & (day_idx <= wet_day_hi), P_W * m_per_h_to_cm_per_d, 0.0)
        pet_d_cmd = np.where(
            (day_idx > wet_day_lo) & (day_idx <= wet_day_hi),
            E_W * m_per_h_to_cm_per_d,
            np.where((day_idx > wet_day_hi) & (day_idx <= dry_day_hi), E_D * m_per_h_to_cm_per_d, 0.0),
        )
        ptt_d_cmd = np.where(
            (day_idx > wet_day_lo) & (day_idx <= wet_day_hi),
            T_W * m_per_h_to_cm_per_d,
            np.where((day_idx > wet_day_hi) & (day_idx <= dry_day_hi), T_D * m_per_h_to_cm_per_d, 0.0),
        )
        summary = dict(
            T_WU_h=T_WU,
            T_WET_h=T_WET,
            T_DRY_h=T_DRY,
            T_COOL_h=T_COOL,
            T_TOTAL_d=T_TOTAL_d,
            P_W_mh=P_W,
            E_W_mh=E_W,
            T_W_mh=T_W,
            E_D_mh=E_D,
            T_D_mh=T_D,
            P_W_cmh=P_W * 100.0,
            E_W_cmh=E_W * 100.0,
            T_W_cmh=T_W * 100.0,
            E_D_cmd=E_D * m_per_h_to_cm_per_d,
            T_D_cmd=T_D * m_per_h_to_cm_per_d,
        )

    # Vegetation ------------------------------------------------------------
    veg = sot["vegetation"][setup]
    summary["root_d_initial"] = veg["root_d_initial"]
    summary["root_d_final"] = veg["root_d_final"]
    summary["soil_layers"] = layers
    summary["Z_TOP"] = Z_TOP
    summary["Z_BOT"] = Z_BOT
    summary["DZ"] = DZ
    summary["NL"] = NL
    summary["Z_WT"] = Z_WT

    # Tuned params ----------------------------------------------------------
    tuned = sot.get("oat_tuned", {}).get(soil_tag, {}).get(setup, {})

    # Output paths ----------------------------------------------------------
    outputs_root = HERE / "outputs" / setup / soil_tag
    hydrus_ws = outputs_root / "phydrus" / "hydrus1d"

    return ExperimentSpec(
        setup=setup,
        soil_tag=soil_tag,
        Z_TOP=Z_TOP,
        Z_BOT=Z_BOT,
        DZ=DZ,
        NL=NL,
        Z_WT=Z_WT,
        BNDS=BNDS,
        LAYER_DEPTH_CM=LAYER_DEPTH_CM,
        materials=materials,
        soil_db=soil_db,
        layer_mat_ids=layer_mat_ids,
        theta_r=theta_r,
        theta_s=theta_s,
        T_TOTAL_d=T_TOTAL_d,
        phase_bounds_d=phase_bounds_d,
        phase_labels=phase_labels,
        phase_colors=phase_colors,
        t_d=t_d,
        n_steps=n_steps,
        dt_unit=dt_unit,
        prec_native=prec_native,
        pet_native=pet_native,
        ptt_native=ptt_native,
        prec_d_cmd=prec_d_cmd,
        pet_d_cmd=pet_d_cmd,
        ptt_d_cmd=ptt_d_cmd,
        day_idx=day_idx,
        root_d_initial=veg["root_d_initial"],
        root_d_final=veg["root_d_final"],
        root_growth_model=veg.get("root_growth", "linear"),
        tuned=tuned,
        outputs_root=outputs_root,
        hydrus_ws=hydrus_ws,
        summary=summary,
    )


# ─────────────────────────────────────────────────────────────────────────────
# HYDRUS-1D reference run
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class HydrusResult:
    t_d: np.ndarray  # plotting axis (matches spec.t_d)
    h_depth_cm: np.ndarray  # bottom-first depth axis [cm]
    gw_depth_cm: np.ndarray  # WT depth, aligned with t_d [cm, +ve below surface]
    sw_cm: np.ndarray  # surface ponding depth, aligned with t_d [cm]
    theta: np.ndarray  # (T, n_nodes), bottom-first, aligned with t_d [m3 m-3]
    cum_ae_cm: np.ndarray
    cum_at_cm: np.ndarray
    mb_err_cm: np.ndarray  # cumulative external mass-balance error [cm]
    wall_s: float


def run_hydrus(
    spec: ExperimentSpec, *, hydrus_exe: Path | None = None, surface_layer: bool | None = None
) -> HydrusResult:
    """Drive HYDRUS-1D for `spec` via phydrus, parse and align outputs.

    surface_layer:
        True  -> top_bc=2 (atmospheric BC with surface layer; ponding accumulates).
        False -> top_bc=3 (surface runoff: excess infiltration shed).
        None  -> default by setup: basic=False (runoff), intensive=True (ponding).
    """
    import phydrus as ps

    if hydrus_exe is None:
        hydrus_exe = HERE / "hydrus1d" / "bin" / "hydrus"
    if surface_layer is None:
        surface_layer = spec.setup == "intensive"

    spec.hydrus_ws.mkdir(parents=True, exist_ok=True)

    Z_COL_CM = (spec.Z_TOP - spec.Z_BOT) * 100.0
    Z_WT_CM = (spec.Z_TOP - spec.Z_WT) * 100.0

    ml = ps.Model(
        exe_name=str(hydrus_exe),
        ws_name=str(spec.hydrus_ws),
        name="gwswex-reference",
        description=f"{spec.T_TOTAL_d}-d {spec.setup}/{spec.soil_tag} column",
        mass_units="-",
        time_unit="days",
        length_unit="cm",
        print_screen=False,
    )
    ml.add_time_info(
        tinit=0,
        tmax=spec.T_TOTAL_d,
        print_times=True,
        dt=0.001 if spec.setup == "intensive" else 0.01,
        dtmin=1e-6,
        dtmax=1 / 24,
        dtprint=1 / 24,
    )
    if surface_layer:
        ml.add_waterflow(top_bc=2, bot_bc=1, rbot=0.0, rroot=0.0, maxit=200, tolh=1.0, tolth=1e-3, ha=1e-6, hb=1e4)
    else:
        ml.add_waterflow(top_bc=3, bot_bc=1, rbot=0.0, rroot=0.0, maxit=20, tolh=1.0, ha=1e-6, hb=1e4)

    # Materials in cm units
    mat_df = ml.get_empty_material_df(n=len(spec.materials))
    for mat in spec.materials:
        s = spec.soil_db[mat["name"]]
        mat_df.loc[mat["id"]] = [
            s["theta_r"],
            s["theta_s"],
            s["alpha"] / 100.0,
            s["n"],
            s["K_sat"] * 100.0,
            s["lam"],
        ]
    ml.add_material(mat_df)

    # 1-cm node profile, hydrostatic IC, layered material assignment
    profile = ps.create_profile(top=0.0, bot=-Z_COL_CM, dx=1.0, mat=spec.materials[0]["id"])
    profile["h"] = -Z_WT_CM - profile["x"].to_numpy(float)
    z_node_depth_cm = -profile["x"].to_numpy(float)
    dz_layer_cm = (spec.Z_TOP - spec.Z_BOT) * 100.0 / spec.NL
    node_layer_idx = np.minimum(np.floor(z_node_depth_cm / dz_layer_cm).astype(int), spec.NL - 1)
    profile["Mat"] = np.array([spec.layer_mat_ids[i] for i in node_layer_idx], dtype=int)

    # Linear root density Beta(z) (intensive uses constant root depth)
    z = profile["x"].to_numpy(float)
    root_depth_cm = spec.root_d_initial * 100.0
    beta = np.where(z >= -root_depth_cm, np.maximum(1.0 + z / root_depth_cm, 0.0), 0.0)
    dz_cm = abs(float(profile["x"].iloc[1] - profile["x"].iloc[0]))
    bint = float(np.trapezoid(beta, dx=dz_cm)) if hasattr(np, "trapezoid") else float(np.trapz(beta, dx=dz_cm))
    if bint > 0:
        beta /= bint
    profile["Beta"] = beta

    ml.add_profile(profile)
    ml.add_obs_nodes(
        [-25.0, -50.0, -100.0, -150.0, -200.0, -250.0]
        if spec.setup == "basic"
        else [-25.0, -50.0, -75.0, -100.0, -120.0, -140.0]
    )

    atm = pd.DataFrame(
        {
            "tAtm": spec.day_idx,
            "Prec": spec.prec_d_cmd,
            "rSoil": spec.pet_d_cmd,
            "rRoot": spec.ptt_d_cmd,
            "hCritA": np.full_like(spec.day_idx, 1e5),
            "rB": np.zeros_like(spec.day_idx),
            "hB": np.zeros_like(spec.day_idx),
            "ht": np.zeros_like(spec.day_idx),
        }
    )
    ml.add_atmospheric_bc(atm)
    ml.add_root_uptake(model=0, poptm=[-25] * len(spec.materials), p0=-10, p2h=-200, p2l=-800, p3=-8000)
    ml.add_root_growth(
        irootin=2,
        irfak=1,
        trmin=0,
        trmed=0,
        trmax=spec.T_TOTAL_d,
        xrmin=spec.root_d_initial * 100.0,
        xrmed=0,
        xrmax=spec.root_d_final * 100.0,
        trperiod=365,
    )
    ml.write_input()

    t0 = time.perf_counter()
    result = ml.simulate()
    wall_s = time.perf_counter() - t0
    if result.returncode != 0:
        raise RuntimeError(f"HYDRUS-1D failed - check {spec.hydrus_ws / 'Error.msg'}")

    # ── Parse NOD_INF ────────────────────────────────────────────────────
    nod_raw = ml.read_nod_inf()
    if isinstance(nod_raw, pd.DataFrame):
        nod_raw = {0.0: nod_raw}

    times, h_list, sm_list, depth_ref = [], [], [], None
    for tk in sorted(nod_raw.keys()):
        df = nod_raw[tk].dropna(subset=["Depth", "Head", "Moisture"]).sort_values("Depth")
        if depth_ref is None:
            depth_ref = np.abs(df["Depth"].to_numpy(float))
        times.append(tk)
        h_list.append(df["Head"].to_numpy(float))
        sm_list.append(df["Moisture"].to_numpy(float))
    h_t = np.array(times)
    h_depth = depth_ref
    h_theta = np.array(sm_list)
    h_head = np.array(h_list)

    h_gw_all = np.array(
        [
            _interp_wt(
                nod_raw[tk].dropna(subset=["Depth", "Head"]).sort_values("Depth")["Depth"].to_numpy(float),
                nod_raw[tk].dropna(subset=["Depth", "Head"]).sort_values("Depth")["Head"].to_numpy(float),
                allow_pond=surface_layer,
            )
            for tk in sorted(nod_raw.keys())
        ]
    )
    if not surface_layer:
        h_gw_all = np.maximum(h_gw_all, 0.0)
    h_sw_all = np.maximum(h_head[:, -1], 0.0)

    # NaN-mask past the last HYDRUS-reported time so that early termination
    # (Picard nonconvergence in persistently-saturated columns under flux-top
    # plus Dirichlet-bottom BC; cf. intensive sand and intensive sand-loam) is
    # surfaced honestly rather than hidden by np.interp's right-edge clamp.
    def _interp_truncated(xp_target, xp_src, fp_src):
        out = np.interp(xp_target, xp_src, fp_src)
        if xp_src.size > 0:
            out = np.where(xp_target > xp_src[-1], np.nan, out)
        return out

    gw_depth = np.maximum(_interp_truncated(spec.t_d, h_t, h_gw_all), 0.0)
    sw_cm = _interp_truncated(spec.t_d, h_t, h_sw_all)
    theta_d = np.column_stack([_interp_truncated(spec.t_d, h_t, h_theta[:, j]) for j in range(h_theta.shape[1])])

    # ── Cumulative ET + MB from T_LEVEL.OUT ─────────────────────────────
    raw_t, cum_e, cum_t, vol, cum_vtop, cum_vbot = [], [], [], [], [], []
    with open(spec.hydrus_ws / "T_LEVEL.OUT") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 22:
                continue
            try:
                raw_t.append(float(parts[0]))
                cum_vtop.append(float(parts[8]))
                cum_t.append(float(parts[9]))
                cum_vbot.append(float(parts[10]))
                vol.append(float(parts[16]))
                cum_e.append(float(parts[18]))
            except ValueError:
                continue
    raw_t = np.array(raw_t)
    cum_e = np.array(cum_e)
    cum_t = np.array(cum_t)
    vol = np.array(vol)
    cum_vtop = np.array(cum_vtop)
    cum_vbot = np.array(cum_vbot)
    cum_ae_d = _interp_truncated(spec.t_d, raw_t, cum_e)
    cum_at_d = _interp_truncated(spec.t_d, raw_t, cum_t)
    vol_d = _interp_truncated(spec.t_d, raw_t, vol)
    vtop_d = _interp_truncated(spec.t_d, raw_t, cum_vtop)
    vbot_d = _interp_truncated(spec.t_d, raw_t, cum_vbot)
    mb_err_d = vol_d - vol[0] + vtop_d - vbot_d + cum_at_d

    return HydrusResult(
        t_d=spec.t_d,
        h_depth_cm=h_depth,
        gw_depth_cm=gw_depth,
        sw_cm=sw_cm,
        theta=theta_d,
        cum_ae_cm=cum_ae_d,
        cum_at_cm=cum_at_d,
        mb_err_cm=mb_err_d,
        wall_s=wall_s,
    )


def _interp_wt(z_asc: np.ndarray, h_asc: np.ndarray, *, allow_pond: bool = False) -> float:
    sc = np.where(np.diff(np.signbit(h_asc)))[0]
    if sc.size == 0 or np.all(h_asc > 0) or np.all(h_asc < 0):
        if allow_pond and np.all(h_asc >= 0):
            return -float(z_asc[-1]) - float(h_asc[-1])
        return np.nan
    i = sc[-1]
    dh = h_asc[i + 1] - h_asc[i]
    z0 = z_asc[i] if abs(dh) < 1e-12 else z_asc[i] + (-h_asc[i]) * (z_asc[i + 1] - z_asc[i]) / dh
    return -z0


# ─────────────────────────────────────────────────────────────────────────────
# GWSWEX result aggregation (used by notebook collect callback)
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class GwswexResult:
    label: str
    gw_depth_cm: np.ndarray  # +ve below surface, aligned with t_d
    sw_cm: np.ndarray
    theta: np.ndarray  # (T, NL) top-to-bottom
    cum_ae_cm: np.ndarray
    cum_at_cm: np.ndarray
    int_mb_err_cm: np.ndarray
    ext_mb_err_cm: np.ndarray | None = None
    wall_s: float = 0.0


def collect_basic(spec: ExperimentSpec, m, *, steps_per_day: int = 1) -> GwswexResult:
    """Run a basic-setup (daily/hourly) model with model.run() callback.

    The callback collects per-step state, mass-balance and intrinsic MB error.
    For implicit (steps_per_day=24) it sub-samples to end-of-day for plotting.
    """
    gw_all, sw_all, th_all, ae_all, at_all, mb_step = [], [], [], [], [], []

    def _cb(_t, state):
        gw_all.append(state["GWH"][0])
        sw_all.append(float(state["SW"][0]))
        th_all.append(state["theta"][:, 0].copy())
        mb = m.get_mass_balance()
        ae_all.append(float(mb["evap"][0]))
        at_all.append(float(mb["transp"][0]))
        ds = float(mb["delta_gw"][0]) + float(mb["delta_sw"][0]) + float(mb["delta_uz"][0])
        net = float(mb["precip"][0]) - float(mb["evap"][0]) - float(mb["transp"][0])
        mb_step.append(ds - net)

    t0 = time.perf_counter()
    m.run(callback=_cb)
    wall = time.perf_counter() - t0
    m.deinit()

    gw_arr = np.array(gw_all)
    sw_arr = np.array(sw_all)
    th_arr = np.array(th_all)
    cum_ae = np.cumsum(ae_all) * 100.0
    cum_at = np.cumsum(at_all) * 100.0
    int_err = np.cumsum(mb_step) * 100.0
    idx = np.arange(steps_per_day - 1, len(gw_arr), steps_per_day)
    return GwswexResult(
        label=getattr(m, "name", "gwswex"),
        gw_depth_cm=(spec.Z_TOP - gw_arr[idx]) * 100.0,
        sw_cm=sw_arr[idx] * 100.0,
        theta=th_arr[idx],
        cum_ae_cm=cum_ae[idx],
        cum_at_cm=cum_at[idx],
        int_mb_err_cm=int_err[idx],
        wall_s=wall,
    )


def collect_intensive(spec: ExperimentSpec, m) -> GwswexResult:
    """Step a 1-h-macro model (intensive setup), applying hourly forcing.

    Records both intrinsic and external (storage-reconstruction) MB error.
    """
    N = spec.n_steps
    gwh = np.zeros(N)
    sw = np.zeros(N)
    uz = np.zeros(N)
    gwv = np.zeros(N)
    theta = np.zeros((N, spec.NL))
    p = np.zeros(N)
    e = np.zeros(N)
    tr = np.zeros(N)
    dgw = np.zeros(N)
    dsw = np.zeros(N)
    duz = np.zeros(N)

    s0 = m.get_state()
    z_sat0 = float(np.clip(float(s0["GWH"][0]) - spec.Z_BOT, 0.0, spec.Z_TOP - spec.Z_BOT))
    ph0 = float(s0["GWV"][0]) + float(np.sum(s0["UZ"][:, 0])) + float(s0["SW"][0]) + spec.theta_r * z_sat0

    t0 = time.perf_counter()
    for t in range(N):
        m.update_forcing(
            t, precip=float(spec.prec_native[t]), pet=float(spec.pet_native[t]), ptt=float(spec.ptt_native[t])
        )
        m.run_step(t)
        st = m.get_state()
        mb = m.get_mass_balance()
        gwh[t] = float(st["GWH"][0])
        sw[t] = float(st["SW"][0])
        uz[t] = float(np.sum(st["UZ"][:, 0]))
        gwv[t] = float(st["GWV"][0])
        theta[t] = st["theta"][:, 0].copy()
        p[t] = float(mb["precip"][0])
        e[t] = float(mb["evap"][0])
        tr[t] = float(mb["transp"][0])
        dgw[t] = float(mb["delta_gw"][0])
        dsw[t] = float(mb["delta_sw"][0])
        duz[t] = float(mb["delta_uz"][0])
    wall = time.perf_counter() - t0
    m.deinit()

    z_sat = np.clip(gwh - spec.Z_BOT, 0.0, spec.Z_TOP - spec.Z_BOT)
    physical = (gwv + uz + sw + spec.theta_r * z_sat) * 100.0
    ext_err = physical - ph0 * 100.0 - np.cumsum(p - e - tr) * 100.0
    int_err = np.cumsum(((dgw + dsw + duz) - (p - e - tr))) * 100.0

    return GwswexResult(
        label=getattr(m, "name", "gwswex"),
        gw_depth_cm=(spec.Z_TOP - gwh) * 100.0,
        sw_cm=sw * 100.0,
        theta=theta,
        cum_ae_cm=np.cumsum(e) * 100.0,
        cum_at_cm=np.cumsum(tr) * 100.0,
        int_mb_err_cm=int_err,
        ext_mb_err_cm=ext_err,
        wall_s=wall,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

# Tol muted palette (CVD-safe). Consistent across all comparison notebooks:
#   HYDRUS-1D = indigo, GWSWEX implicit = wine, GWSWEX explicit = forest green.
STYLE = {
    "HYDRUS-1D": dict(color="#332288", lw=2.0, ls="-", label="HYDRUS-1D"),
    "GWSWEX implicit": dict(color="#882255", lw=1.8, ls="-", label="GWSWEX implicit"),
    "GWSWEX explicit": dict(color="#117733", lw=1.8, ls="--", label="GWSWEX explicit"),
}


def setup_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "figure.dpi": 130,
        }
    )


def shade_phases(ax, spec: ExperimentSpec, *, label: bool = True) -> None:
    ylo, yhi = ax.get_ylim()
    for x0, x1, lbl, col in zip(
        spec.phase_bounds_d[:-1], spec.phase_bounds_d[1:], spec.phase_labels, spec.phase_colors
    ):
        ax.axvspan(x0, x1, color=col, alpha=0.55, zorder=0)
        if label:
            ax.text(0.5 * (x0 + x1), yhi - 0.04 * (yhi - ylo), lbl, fontsize=8, ha="center", va="top", alpha=0.75)
    ax.set_ylim(ylo, yhi)


# ─── Metrics ────────────────────────────────────────────────────────────────


def rmse(sim, obs):
    sim = np.asarray(sim, float)
    obs = np.asarray(obs, float)
    m = np.isfinite(sim) & np.isfinite(obs)
    if m.sum() == 0:
        return float("nan")
    return float(np.sqrt(np.mean((sim[m] - obs[m]) ** 2)))


def mae(sim, obs):
    sim = np.asarray(sim, float)
    obs = np.asarray(obs, float)
    m = np.isfinite(sim) & np.isfinite(obs)
    if m.sum() == 0:
        return float("nan")
    return float(np.mean(np.abs(sim[m] - obs[m])))


def bias(sim, obs):
    sim = np.asarray(sim, float)
    obs = np.asarray(obs, float)
    m = np.isfinite(sim) & np.isfinite(obs)
    if m.sum() == 0:
        return float("nan")
    return float(np.mean(sim[m] - obs[m]))


def nse(sim, obs):
    sim = np.asarray(sim, float)
    obs = np.asarray(obs, float)
    m = np.isfinite(sim) & np.isfinite(obs)
    if m.sum() == 0:
        return float("nan")
    o = obs[m]
    s = sim[m]
    ss = np.sum((o - o.mean()) ** 2)
    return float(1.0 - np.sum((s - o) ** 2) / ss) if ss > 0 else float("nan")


# ─── Standard plot panels ───────────────────────────────────────────────────


def plot_water_table(spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult]) -> None:
    """3-panel: WT depth (top, big), ponded SW (middle), residual (bottom)."""
    fig, (ax_wt, ax_sw, ax_res) = plt.subplots(
        3,
        1,
        figsize=(10, 7.0),
        gridspec_kw={"height_ratios": [3, 1, 1]},
        sharex=True,
    )
    fig.suptitle(f"Water-table depth: {spec.T_TOTAL_d}-d {spec.soil_tag} column " f"({spec.setup})", fontsize=11)

    Z_WT_CM = (spec.Z_TOP - spec.Z_WT) * 100.0
    ax_wt.plot(spec.t_d, hyd.gw_depth_cm, **STYLE["HYDRUS-1D"])
    for key, res in runs.items():
        ax_wt.plot(spec.t_d, np.maximum(res.gw_depth_cm, 0.0), **STYLE[key])
    ax_wt.axhline(Z_WT_CM, color="k", ls="--", lw=0.8, alpha=0.4, label=f"Initial WT ({Z_WT_CM:.0f} cm)")
    ax_wt.invert_yaxis()
    ax_wt.set_ylabel("WT depth [cm]")
    ax_wt.legend(loc="lower right", fontsize=9)
    shade_phases(ax_wt, spec)

    ax_sw.plot(spec.t_d, hyd.sw_cm, **STYLE["HYDRUS-1D"])
    for key, res in runs.items():
        ax_sw.plot(spec.t_d, res.sw_cm, lw=1.4, color=STYLE[key]["color"], ls=STYLE[key]["ls"], label=key)
    sw_max = max(
        float(np.nanmax(hyd.sw_cm)) if np.any(np.isfinite(hyd.sw_cm)) else 0.0,
        *(float(np.nanmax(r.sw_cm)) if np.any(np.isfinite(r.sw_cm)) else 0.0 for r in runs.values()),
        1e-3,
    )
    ax_sw.set_ylabel("Ponded SW [cm]")
    ax_sw.axhline(0, color="k", lw=0.5, alpha=0.5)
    ax_sw.set_ylim(-0.05 * sw_max, 1.1 * sw_max)
    ax_sw.legend(loc="upper right", fontsize=8)
    shade_phases(ax_sw, spec, label=False)

    for key, res in runs.items():
        resid = hyd.gw_depth_cm - np.maximum(res.gw_depth_cm, 0.0)
        ax_res.fill_between(spec.t_d, resid, 0, alpha=0.25, color=STYLE[key]["color"])
        ax_res.plot(spec.t_d, resid, lw=1.2, color=STYLE[key]["color"], label=key)
    ax_res.axhline(0, color="k", lw=0.5)
    ax_res.set_ylabel("DIFF [cm]")
    ax_res.set_xlabel("Time [days]")
    ax_res.legend(loc="lower left", fontsize=8)
    shade_phases(ax_res, spec, label=False)
    fig.tight_layout()
    plt.show()

    for key, res in runs.items():
        sim = np.maximum(res.gw_depth_cm, 0.0)
        print(
            f"{key}: RMSE={rmse(sim, hyd.gw_depth_cm):.2f} cm, "
            f"MAE={mae(sim, hyd.gw_depth_cm):.2f} cm, "
            f"NSE={nse(sim, hyd.gw_depth_cm):.3f}, "
            f"Bias={bias(sim, hyd.gw_depth_cm):+.2f} cm"
        )


def plot_zone_theta(
    spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult], zones: list[tuple[float, float, str]]
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True)
    fig.suptitle("Zone-averaged volumetric water content", fontsize=11)
    for ax, (z1, z2, label) in zip(axes.flat, zones):
        hm = (hyd.h_depth_cm >= z1) & (hyd.h_depth_cm <= z2)
        h_zone = hyd.theta[:, hm].mean(axis=1)
        gm = (spec.LAYER_DEPTH_CM >= z1) & (spec.LAYER_DEPTH_CM <= z2)
        ax.plot(spec.t_d, h_zone, **STYLE["HYDRUS-1D"])
        for key, res in runs.items():
            ax.plot(spec.t_d, res.theta[:, gm].mean(1), **STYLE[key])
        vals = np.concatenate([h_zone, *(r.theta[:, gm].mean(1) for r in runs.values())])
        ymin, ymax = float(np.nanmin(vals)), float(np.nanmax(vals))
        ax.set_ylim(max(0.0, ymin - 0.05), min(1.0, ymax + 0.05))
        ax.set_title(label, fontsize=9)
        shade_phases(ax, spec, label=False)
        for j, (key, res) in enumerate(runs.items()):
            sim = res.theta[:, gm].mean(1)
            ax.text(
                0.97,
                0.06 + 0.12 * j,
                f"{key.split()[-1][:4]} RMSE={rmse(sim, h_zone):.4f}",
                transform=ax.transAxes,
                ha="right",
                fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", alpha=0.8),
            )
    axes[0, 0].set_ylabel(r"$\theta$ [m$^3$ m$^{-3}$]")
    axes[1, 0].set_ylabel(r"$\theta$ [m$^3$ m$^{-3}$]")
    axes[1, 0].set_xlabel("Time [days]")
    axes[1, 1].set_xlabel("Time [days]")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    plt.show()


def plot_theta_snapshots(
    spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult], snaps: list[tuple[int, str]]
) -> None:
    n = len(snaps)
    fig, axes = plt.subplots(1, n, figsize=(2.5 * n + 2, 5.5), sharey=True)
    if n == 1:
        axes = [axes]
    fig.suptitle(r"$\theta$ profiles at phase snapshots", fontsize=11)
    for ax, (idx, label) in zip(axes, snaps):
        ax.plot(hyd.theta[idx, :], hyd.h_depth_cm, **STYLE["HYDRUS-1D"])
        ax.axhline(hyd.gw_depth_cm[idx], color=STYLE["HYDRUS-1D"]["color"], ls="--", lw=0.9, alpha=0.5)
        for key, res in runs.items():
            ax.plot(res.theta[idx, :], spec.LAYER_DEPTH_CM, **STYLE[key])
            ax.axhline(max(res.gw_depth_cm[idx], 0.0), color=STYLE[key]["color"], ls="--", lw=0.9, alpha=0.5)
        ax.set_title(label, fontsize=9)
        ax.set_xlabel(r"$\theta$ [m$^3$ m$^{-3}$]")
        ax.set_xlim(0.0, 0.45)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Depth below surface [cm]")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    plt.show()


def plot_cum_et(spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult]) -> None:
    et_labels = ["Actual evaporation (AE)", "Actual transpiration (AT)", "Total AET"]
    h_et = [hyd.cum_ae_cm, hyd.cum_at_cm, hyd.cum_ae_cm + hyd.cum_at_cm]
    sim_et = {k: [r.cum_ae_cm, r.cum_at_cm, r.cum_ae_cm + r.cum_at_cm] for k, r in runs.items()}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))
    fig.suptitle("Cumulative actual evapotranspiration", fontsize=11)
    for ax, lbl, hd in zip(axes, et_labels, h_et):
        ax.plot(spec.t_d, hd, **STYLE["HYDRUS-1D"])
        for key in runs:
            ax.plot(spec.t_d, sim_et[key][et_labels.index(lbl)], **STYLE[key])
        ax.set_title(lbl, fontsize=9)
        ax.set_xlabel("Time [days]")
        ax.set_ylabel("Cumulative [cm]")
        shade_phases(ax, spec, label=False)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    plt.show()

    cols = ["HYDRUS-1D"] + list(runs.keys())
    print(f"{'':>22s}  " + "  ".join(f"{c:>14s}" for c in cols))
    for lbl, hd in zip(["Cum AE [cm]", "Cum AT [cm]", "Cum AET [cm]"], h_et):
        vals = [hd[-1]] + [sim_et[k][["Cum AE [cm]", "Cum AT [cm]", "Cum AET [cm]"].index(lbl)][-1] for k in runs]
        print(f"{lbl:>22s}  " + "  ".join(f"{v:>14.3f}" for v in vals))


def plot_surface_ponding(spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult]) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.0))
    fig.suptitle("Surface ponding depth", fontsize=11)
    ax.plot(spec.t_d, hyd.sw_cm, **STYLE["HYDRUS-1D"])
    for key, res in runs.items():
        ax.plot(spec.t_d, res.sw_cm, **STYLE[key])
    ax.set_ylabel("SW depth [cm]")
    ax.set_xlabel("Time [days]")
    ax.legend(fontsize=8, loc="upper right")
    shade_phases(ax, spec)
    fig.tight_layout()
    plt.show()
    print(f"{'Model':<20s} {'Max SW [cm]':>13s}  {'End SW [cm]':>13s}")
    print("-" * 50)
    print(f"{'HYDRUS-1D':<20s} {hyd.sw_cm.max():>13.2f}  {hyd.sw_cm[-1]:>13.2f}")
    for key, res in runs.items():
        print(f"{key:<20s} {res.sw_cm.max():>13.2f}  {res.sw_cm[-1]:>13.2f}")


def plot_mass_balance(spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult]) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.5))
    fig.suptitle("Cumulative mass-balance error " "(positive = storage gain exceeds boundary input)", fontsize=11)
    ax.plot(spec.t_d, hyd.mb_err_cm, color=STYLE["HYDRUS-1D"]["color"], lw=2.0, ls="-", label="HYDRUS-1D")
    series = [hyd.mb_err_cm]
    for key, res in runs.items():
        ax.plot(spec.t_d, res.int_mb_err_cm, color=STYLE[key]["color"], lw=1.8, ls="-", label=f"{key} intrinsic")
        series.append(res.int_mb_err_cm)
        if res.ext_mb_err_cm is not None:
            ax.plot(spec.t_d, res.ext_mb_err_cm, color=STYLE[key]["color"], lw=1.6, ls="--", label=f"{key} external")
            series.append(res.ext_mb_err_cm)
    ax.axhline(0, color="k", lw=0.6, alpha=0.6)
    ymax = max(1e-3, 1.1 * np.max(np.abs(np.concatenate(series))))
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel("Cumulative error [cm]")
    ax.set_xlabel("Time [days]")
    ax.legend(fontsize=8, loc="best", ncol=2)
    shade_phases(ax, spec, label=False)
    fig.tight_layout()
    plt.show()

    print(f"\nCumulative MB error at t = {spec.T_TOTAL_d} d")
    print(f"{'Model':<24s} {'Final [cm]':>14s}  {'Max |err| [cm]':>16s}")
    print("-" * 60)
    print(f"{'HYDRUS-1D':<24s} {hyd.mb_err_cm[-1]:>+14.6f}  " f"{np.max(np.abs(hyd.mb_err_cm)):>16.6f}")
    for key, res in runs.items():
        print(
            f"{key + ' intrinsic':<24s} {res.int_mb_err_cm[-1]:>+14.6f}  " f"{np.max(np.abs(res.int_mb_err_cm)):>16.6f}"
        )
        if res.ext_mb_err_cm is not None:
            print(
                f"{key + ' external':<24s} {res.ext_mb_err_cm[-1]:>+14.6f}  "
                f"{np.max(np.abs(res.ext_mb_err_cm)):>16.6f}"
            )


def plot_metrics_summary(
    spec: ExperimentSpec, hyd: HydrusResult, runs: dict[str, GwswexResult], zone_specs: list[tuple[float, float, str]]
) -> None:
    """Grouped bar charts: WT-depth metrics by phase + zone-θ metrics by zone."""
    # Phase masks
    pm = []
    pm.append((f"Full\n(0-{spec.T_TOTAL_d} d)", np.ones_like(spec.t_d, dtype=bool)))
    for x0, x1, lbl in zip(spec.phase_bounds_d[:-1], spec.phase_bounds_d[1:], spec.phase_labels):
        pm.append(
            (f"{lbl}\n({int(x0)}-{int(x1)} d)", (spec.t_d > x0) & (spec.t_d <= x1) if x0 > 0 else (spec.t_d <= x1))
        )

    wt_metrics = {}
    for key, res in runs.items():
        wt_metrics[key] = {
            p: dict(
                RMSE=rmse(np.maximum(res.gw_depth_cm, 0.0)[mask], hyd.gw_depth_cm[mask]),
                MAE=mae(np.maximum(res.gw_depth_cm, 0.0)[mask], hyd.gw_depth_cm[mask]),
                NSE=nse(np.maximum(res.gw_depth_cm, 0.0)[mask], hyd.gw_depth_cm[mask]),
                Bias=bias(np.maximum(res.gw_depth_cm, 0.0)[mask], hyd.gw_depth_cm[mask]),
            )
            for p, mask in pm
        }

    th_metrics = {}
    for key, res in runs.items():
        th_metrics[key] = {}
        for z1, z2, zlbl in zone_specs:
            hm = (hyd.h_depth_cm >= z1) & (hyd.h_depth_cm <= z2)
            gm = (spec.LAYER_DEPTH_CM >= z1) & (spec.LAYER_DEPTH_CM <= z2)
            h_zone = hyd.theta[:, hm].mean(axis=1)
            sim = res.theta[:, gm].mean(axis=1)
            th_metrics[key][zlbl] = dict(
                RMSE=rmse(sim, h_zone),
                MAE=mae(sim, h_zone),
                NSE=nse(sim, h_zone),
                Bias=bias(sim, h_zone),
            )

    solver_colors = {k: STYLE[k]["color"] for k in runs}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle("Performance metrics vs HYDRUS-1D", fontsize=12)
    plabels = [p[0] for p in pm]
    x = np.arange(len(plabels))
    w = 0.35
    for col, metric in enumerate(["RMSE", "MAE", "NSE"]):
        ax = axes[0, col]
        for i, (slbl, sc) in enumerate(solver_colors.items()):
            vals = [wt_metrics[slbl][p][metric] for p in plabels]
            bars = ax.bar(x + (i - 0.5) * w, vals, w, color=sc, alpha=0.8, label=slbl)
            for bar, v in zip(bars, vals):
                fmt = f"{v:.3f}" if metric == "NSE" else f"{v:.2f}"
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), fmt, ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels(plabels, fontsize=7.5)
        unit = "" if metric == "NSE" else " [cm]"
        ax.set_title(f"WT depth {metric}{unit}", fontsize=10)
        if col == 0:
            ax.legend(fontsize=8)
    zlabels = [z[2] for z in zone_specs]
    x2 = np.arange(len(zlabels))
    for col, metric in enumerate(["RMSE", "MAE", "NSE"]):
        ax = axes[1, col]
        for i, (slbl, sc) in enumerate(solver_colors.items()):
            vals = [th_metrics[slbl][z][metric] for z in zlabels]
            bars = ax.bar(x2 + (i - 0.5) * w, vals, w, color=sc, alpha=0.8, label=slbl)
            for bar, v in zip(bars, vals):
                fmt = f"{v:.3f}" if metric == "NSE" else f"{v:.5f}"
                ax.text(
                    bar.get_x() + bar.get_width() / 2, bar.get_height(), fmt, ha="center", va="bottom", fontsize=6.5
                )
        ax.set_xticks(x2)
        ax.set_xticklabels(zlabels, fontsize=7.5)
        unit = "" if metric == "NSE" else r" [m$^3$m$^{-3}$]"
        ax.set_title(f"Zone-avg θ {metric}{unit}", fontsize=10)
        ax.set_xlabel("Depth zone")
    fig.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Markdown header rendering (so README values track the JSON SoT)
# ─────────────────────────────────────────────────────────────────────────────


def markdown_header(spec: ExperimentSpec) -> str:
    """Render a Markdown header that mirrors values from experiment_definitions.json.

    Returned string is the body of a notebook markdown cell.
    """
    s = spec.summary
    soil_layers_str = ", ".join(f"{name} ({frac * 100:.0f}%)" for name, frac in s["soil_layers"])
    z_wt_cm = (spec.Z_TOP - spec.Z_WT) * 100.0

    if spec.setup == "basic":
        title = f"# GWSWEX vs HYDRUS-1D - basic / {spec.soil_tag}"
        body = f"""{title}

{spec.T_TOTAL_d}-day synthetic column experiment comparing two GWSWEX solvers
against HYDRUS-1D. All values shown here are loaded directly from
`experiment_definitions.json` at the start of each cell - this header is
hand-mirrored for readability and is verified by the test harness to track
the JSON source of truth.

| Parameter | Value |
|-----------|-------|
| Column | {spec.Z_TOP - spec.Z_BOT:.1f} m, {spec.NL} layers @ {spec.DZ * 100:.1f} cm |
| Soil profile (top→bottom) | {soil_layers_str} |
| Initial WT depth | {z_wt_cm:.0f} cm below surface |
| Warmup (days {s["T_P1"]}-{s["T_P2"]}) | zero forcing |
| Wet (days {s["T_P2"]}-{s["T_P3"]}) | P = {s["P_WET_cmd"]:g}, PE = {s["PE_WET_cmd"]:g}, PT = {s["PT_WET_cmd"]:g} cm d⁻¹ |
| Dry (days {s["T_P3"]}-{s["T_TOTAL"]}) | P = 0, PE = {s["PE_DRY_cmd"]:g}, PT = {s["PT_DRY_cmd"]:g} cm d⁻¹ |
| Root growth | linear, {s["root_d_initial"]:g} → {s["root_d_final"]:g} m |

**GWSWEX explicit**: operator-split cascade with CFL sub-stepping (daily macro-steps).
**GWSWEX implicit**: mixed-form Richards + Picard/TDMA (hourly macro-steps).

Solver- and material-specific tuned parameters (`MODEL_PARAMS`, `ET_STRESS`,
`solver kwargs`) are taken from the `oat_tuned` section of the same JSON file.
Outputs (HYDRUS-1D workspace, NetCDF state files) land under
`outputs/{spec.setup}/{spec.soil_tag}/`.
"""
    else:  # intensive
        title = f"# GWSWEX vs HYDRUS-1D - intensive / {spec.soil_tag}"
        body = f"""{title}

{spec.T_TOTAL_d}-day, 4-phase column experiment with hourly macro-steps and
explicit ponding resolution (HYDRUS-1D top BC = atmospheric-with-surface-layer,
GWSWEX SW reservoir). All values below come from `experiment_definitions.json`.

| Parameter | Value |
|-----------|-------|
| Column | {spec.Z_TOP - spec.Z_BOT:.1f} m, {spec.NL} layers @ {spec.DZ * 100:.1f} cm |
| Soil profile (top→bottom) | {soil_layers_str} |
| Initial WT depth | {z_wt_cm:.0f} cm below surface |
| Warmup ({s["T_WU_h"]} h) | zero forcing |
| Wet ({s["T_WET_h"]} h) | P = {s["P_W_cmh"]:.3f}, PE = {s["E_W_cmh"]:.3f}, PT = {s["T_W_cmh"]:.3f} cm h⁻¹ |
| Dry ({s["T_DRY_h"]} h) | P = 0, PE = {s["E_D_cmd"]:.2f}, PT = {s["T_D_cmd"]:.2f} cm d⁻¹ |
| Cooldown ({s["T_COOL_h"]} h) | zero forcing |
| Root depth (constant) | {s["root_d_initial"]:g} m |

Solver- and material-specific tuned parameters (`MODEL_PARAMS`, `ET_STRESS`,
`solver kwargs`) are taken from the `oat_tuned` section of the same JSON file.
Outputs land under `outputs/{spec.setup}/{spec.soil_tag}/`.
"""
    return body

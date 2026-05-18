"""GWSWEX demo – deep-GWT → ponding → shallow-GWT column experiment (IMPLICIT solver).

Single loam column (ne=1, nl=15, dz=0.1 m):

  Phase     Duration  Forcing
  --------  --------  -------------------------------------------------------
  Warmup     3 d      Zero forcing; hydrostatic equilibrium.
  Wet       10 d      P = 1.25 mm h⁻¹; GWT rises and exceeds the surface.
  Dry       12 d      ET only (0.4 + 0.3 mm h⁻¹); GWT declines to shallow.
  Cooldown   7 d      Zero forcing; free drainage to hydrostatic equilibrium.

Forcing injected phase-by-phase via model.update_forcing() in a step loop.
Solver: implicit (Picard iteration, picard_tol=1e-5, picard_max_iter=100, dt = 1 h).

Mass-balance sketch
  Sy = 0.352. Hydrostatic UZ ~ 0.363 m; max = 0.516 m; deficit ~ 0.153 m.
  Net wet input = (1.25-0.3) mm/h = 0.95 mm/h.
  GWT hits surface in ~153/0.95 ~ 161 h (day ~6.7); SW builds to ~400 mm.
  Dry ET = 0.7 mm/h x 288 h = 202 mm -> GWT drops ~0.57 m below surface.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import numpy as np

from gwswex import GWSWEXmodel

# ── Output directory ──────────────────────────────────────────────────────────
HERE = Path(__file__).parent
OUTDIR = HERE
# ── Debug mode ────────────────────────────────────────────────────────────────
_cli = argparse.ArgumentParser(add_help=False)
_cli.add_argument("--debug", action="store_true", help="Write per-step diagnostics to OUTDIR/debug.csv")
DEBUG = _cli.parse_known_args()[0].debug
# ── Soil: loam ────────────────────────────────────────────────────────────────
THETA_R, THETA_S = 0.078, 0.43
K_SAT = 0.2496 / 24  # 0.0104 m h⁻¹  (from 0.2496 m d⁻¹)
VAN_G = dict(alpha=3.6, n=1.56, theta_r=THETA_R, theta_s=THETA_S)
LAM = 0.5

# ── Column geometry: 1.5 m total, 15 x 0.1 m layers ─────────────────────────
Z_TOP, Z_BOT, NL, DZ, NE = 1.5, 0.0, 15, 0.1, 1
bot = [list(np.round(np.linspace(Z_TOP - DZ, Z_BOT, NL), 10))]

# Layer boundary depths below surface [m]: [0.0, 0.1, ..., 1.5]
DEPTH_BNDS = np.linspace(0.0, Z_TOP, NL + 1)
DEPTH_MIDS = 0.5 * (DEPTH_BNDS[:-1] + DEPTH_BNDS[1:])

# ── Phase lengths [h] ─────────────────────────────────────────────────────────
T_WU, T_WET, T_DRY, T_COOL = 72, 240, 288, 168  # T_COOL: 7 d = 168 h
N = T_WU + T_WET + T_DRY + T_COOL  # 768 h total (32 d)

# ── Atmospheric forcing [m h⁻¹] ───────────────────────────────────────────────
# Wet phase: native (uncompensated) forcings.
P_W, E_W, T_W = 1.25e-3, 0.2e-3, 0.1e-3  # wet: P, PET, PTT [m h^-1]
E_D, T_D = 0.4e-3, 0.3e-3  # dry: PET, PTT

# ── Build model ───────────────────────────────────────────────────────────────
nc_path = OUTDIR / "demo-implicit.nc"
m = GWSWEXmodel(name="demo-implicit", T="h", L="m", output_fpath=nc_path)
m.init_space(ne=NE, nl=NL, top=[[Z_TOP]], bot=bot, sID=[[1] * NL], vID=[[1]])
m.add_material(id=1, name="loam", K_sat=K_SAT, lam=LAM, vanG=VAN_G)
m.add_vegetation(
    id=1,
    name="pasture",
    root_depth_initial=0.6,
    root_depth_final=0.6,
    root_growth_model="linear",
    et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
)

m.init_time(n_steps=N, dt=1.0, dt_min=1 / 60)

m.set_initial_conditions(
    gw=[0.3],  # GWT at 0.3 m elevation -> 1.2 m below surface
    sw=[0.0],
    uz=[-999],  # hydrostatic equilibrium
)

m.set_solver(solver="implicit", picard_tol=1e-5, picard_max_iter=100)
m.set_model_params(psi_f=0.09, F_min=1e-6, ICratio_min=0.05)

m.set_forcing(precip=0.0, pet=0.0, ptt=0.0)
m.init()  # opens the built-in NetCDF writer at nc_path

# ── Snapshots ────────────────────────────────────────────────────────────────
SNAPS = {
    "After warmup": T_WU - 1,
    "Mid-wet": T_WU + T_WET // 2 - 1,
    "End wet / start dry": T_WU + T_WET - 1,
    "Mid-dry": T_WU + T_WET + T_DRY // 2 - 1,
    "End dry": T_WU + T_WET + T_DRY - 1,
    "End cooldown": N - 1,
}
SNAP_IDX = set(SNAPS.values())
SNAP_ORDER = ["Init (t=0)"] + list(SNAPS.keys())

# ── Time-series storage ───────────────────────────────────────────────────────
ts_gwh = np.zeros(N)
ts_sw = np.zeros(N)
ts_uz = np.zeros(N)
ts_sat = np.zeros(N)  # column-avg relative saturation S=(theta-theta_r)/(theta_s-theta_r)
ts_mbi = np.zeros(N)  # per-step intrinsic MB imbalance [m]
# Prescribed (potential) and actual fluxes [m h⁻¹] per step
ts_p = np.zeros(N)  # prescribed precip
ts_pe = np.zeros(N)  # prescribed PET
ts_pt = np.zeros(N)  # prescribed PTT
ts_ae = np.zeros(N)  # actual evap (from MB)
ts_at = np.zeros(N)  # actual transp (from MB)

# Profile storage at snapshot times
snap_theta: dict[str, np.ndarray] = {}
snap_gwh: dict[str, float] = {}
snap_sw: dict[str, float] = {}

# ── Initial state ─────────────────────────────────────────────────────────────
s0 = m.get_state()
snap_gwh["Init (t=0)"] = float(s0["GWH"][0])
snap_sw["Init (t=0)"] = float(s0["SW"][0])
# Kernel populates theta_curr at the end of kernel_set_ic (hydrostatic VG profile),
# so we can read the initial moisture profile directly from the model state.
snap_theta["Init (t=0)"] = s0["theta"][:, 0].copy()

print(f"\nColumn:  top = {Z_TOP:.2f} m  |  bottom = {Z_BOT:.2f} m  |  {NL} layers x {DZ:.2f} m")
print(f"Init GWT: {snap_gwh['Init (t=0)']:.3f} m elevation  ({Z_TOP - snap_gwh['Init (t=0)']:.2f} m below surface)\n")

HDR = f"{'Time point':<26s}  {'GWH [m]':>8s}  {'UZ [m]':>8s}  {'SW [m]':>8s}"
SEP = "-" * len(HDR)
print(f"{HDR}\n{SEP}")
print(f"{'Init  (t=0)':26s}  {s0['GWH'][0]:8.4f}  {np.sum(s0['UZ'][:, 0]):8.6f}  {s0['SW'][0]:8.6f}")
# ── Debug output setup ─────────────────────────────────────────────────────────────
_DBG_FIELDS = [
    "step_h",
    "phase",
    "gwh_m",
    "sw_m",
    "uz_m",
    "theta_mean",
    "sat_mean",
    "precip_mh",
    "pet_mh",
    "ptt_mh",
    "precip_acc_m",
    "evap_acc_m",
    "transp_acc_m",
    "runoff_acc_m",
    "infiltration_acc_m",
    "recharge_acc_m",
    "lat_gw_m",
    "lat_sw_m",
    "delta_gw_m",
    "delta_sw_m",
    "delta_uz_m",
    "mbi_m",
    "n_substeps",
]
if DEBUG:
    _dbg_fh = open(OUTDIR / "debug.csv", "w", newline="")
    _dbg_writer = csv.DictWriter(_dbg_fh, fieldnames=_DBG_FIELDS)
    _dbg_writer.writeheader()
    print(f"\n[DEBUG] Writing per-step diagnostics -> {OUTDIR / 'debug.csv'}")
    print(
        f"[DEBUG] {'step':>4s}  {'phase':7s}  {'gwh[m]':>8s}  "
        f"{'sw[mm]':>8s}  {'uz[mm]':>7s}  {'sat':>5s}  {'mbi[m]':>10s}  {'nsub':>4s}"
    )
# ── Step loop ────────────────────────────────────────────────────────────────
for t in m.Time.steps:
    if t < T_WU:
        p_f, pe_f, pt_f = 0.0, 0.0, 0.0
    elif t < T_WU + T_WET:
        p_f, pe_f, pt_f = P_W, E_W, T_W
    elif t < T_WU + T_WET + T_DRY:
        p_f, pe_f, pt_f = 0.0, E_D, T_D
    else:
        p_f, pe_f, pt_f = 0.0, 0.0, 0.0  # cooldown: no forcing

    m.update_forcing(t, precip=p_f, pet=pe_f, ptt=pt_f)
    m.run_step(t)

    s = m.get_state()
    mb = m.get_mass_balance()

    ts_gwh[t] = s["GWH"][0]
    ts_sw[t] = s["SW"][0]
    ts_uz[t] = np.sum(s["UZ"][:, 0])
    ts_sat[t] = float(np.mean((s["theta"][:, 0] - THETA_R) / (THETA_S - THETA_R)))
    ts_p[t] = p_f
    ts_pe[t] = pe_f
    ts_pt[t] = pt_f
    ts_ae[t] = float(mb["evap"][0])  # actual evap this step [m SI = m/step]
    ts_at[t] = float(mb["transp"][0])  # actual transp this step [m SI]

    ds = mb["delta_gw"][0] + mb["delta_sw"][0] + mb["delta_uz"][0]
    # Standalone kernel (no SW routing): rejected precipitation is captured by
    # SW_curr and therefore already appears in delta_sw.  acc_runoff is the
    # depth an external SW router would have to accept — it is NOT a physical
    # outflow from the element, so it must not be subtracted from net.
    net = mb["precip"][0] - mb["evap"][0] - mb["transp"][0]
    ts_mbi[t] = ds - net

    if DEBUG:
        _phase = (
            "warmup" if t < T_WU else "wet" if t < T_WU + T_WET else "dry" if t < T_WU + T_WET + T_DRY else "cooldown"
        )
        print(
            f"[DEBUG] {t+1:4d}  {_phase:7s}  {ts_gwh[t]:8.4f}  {ts_sw[t]*1e3:8.3f}  "
            f"{ts_uz[t]*1e3:7.3f}  {ts_sat[t]:5.3f}  {ts_mbi[t]:10.3e}  {int(mb['n_substeps'][0]):4d}"
        )
        _dbg_writer.writerow(
            {
                "step_h": t + 1,
                "phase": _phase,
                "gwh_m": f"{ts_gwh[t]:.6f}",
                "sw_m": f"{ts_sw[t]:.6e}",
                "uz_m": f"{ts_uz[t]:.6e}",
                "theta_mean": f"{float(np.mean(s['theta'][:, 0])):.6f}",
                "sat_mean": f"{ts_sat[t]:.6f}",
                "precip_mh": f"{p_f:.4e}",
                "pet_mh": f"{pe_f:.4e}",
                "ptt_mh": f"{pt_f:.4e}",
                "precip_acc_m": f"{float(mb['precip'][0]):.6e}",
                "evap_acc_m": f"{float(mb['evap'][0]):.6e}",
                "transp_acc_m": f"{float(mb['transp'][0]):.6e}",
                "runoff_acc_m": f"{float(mb['runoff'][0]):.6e}",
                "infiltration_acc_m": f"{float(mb['infiltration'][0]):.6e}",
                "recharge_acc_m": f"{float(mb['recharge'][0]):.6e}",
                "lat_gw_m": f"{float(mb['lat_gw'][0]):.6e}",
                "lat_sw_m": f"{float(mb['lat_sw'][0]):.6e}",
                "delta_gw_m": f"{float(mb['delta_gw'][0]):.6e}",
                "delta_sw_m": f"{float(mb['delta_sw'][0]):.6e}",
                "delta_uz_m": f"{float(mb['delta_uz'][0]):.6e}",
                "mbi_m": f"{ts_mbi[t]:.6e}",
                "n_substeps": int(mb["n_substeps"][0]),
            }
        )

    if t in SNAP_IDX:
        lbl = next(k for k, v in SNAPS.items() if v == t)
        snap_theta[lbl] = s["theta"][:, 0].copy()
        snap_gwh[lbl] = float(ts_gwh[t])
        snap_sw[lbl] = float(ts_sw[t])
        tag = f"t={t + 1:>4d} h  {lbl}"
        print(f"{tag:<26s}  {ts_gwh[t]:8.4f}  {ts_uz[t]:8.6f}  {ts_sw[t]:8.6f}")

print(f"Saved: {nc_path}")
m.deinit()  # closes the built-in NetCDF writer
print(SEP)

if DEBUG:
    _dbg_fh.close()
    print(f"[DEBUG] Saved: {OUTDIR / 'debug.csv'}")

# ── Derived cumulative series [mm] ────────────────────────────────────────────
# mb fluxes are in SI (m per step); forcing rates are in m h⁻¹ * 1 h/step = m/step
t_d = np.arange(1, N + 1) / 24.0
d_wu = T_WU / 24
d_dry = (T_WU + T_WET) / 24
d_cool = (T_WU + T_WET + T_DRY) / 24
mb_cum = np.cumsum(ts_mbi) * 1e3  # [mm]
cum_p = np.cumsum(ts_p) * 1e3  # prescribed P [mm]
cum_pe = np.cumsum(ts_pe) * 1e3  # prescribed PE [mm]
cum_pt = np.cumsum(ts_pt) * 1e3  # prescribed PT [mm]
# MB fluxes are in SI (m s⁻¹ * s = m); convert to mm.
# The model uses T="h", so mb values are in m per step (1 h).
# They are stored internally in SI (m s⁻¹); get_mass_balance returns SI.
# One step = 1 h = 3600 s; mb["evap"] is in m (volume per unit area per step in SI).
cum_ae = np.cumsum(ts_ae) * 1e3  # actual AE [mm]
cum_at = np.cumsum(ts_at) * 1e3  # actual AT [mm]

# ── Shared plot helpers ───────────────────────────────────────────────────────
PCOL = {"Warmup": "#ddeef7", "Wet": "#ceecc8", "Dry": "#faecd1", "Cooldown": "#e8e0f0"}


def _shade(ax):
    ax.axvspan(0, d_wu, color=PCOL["Warmup"], alpha=0.85, zorder=0)
    ax.axvspan(d_wu, d_dry, color=PCOL["Wet"], alpha=0.85, zorder=0)
    ax.axvspan(d_dry, d_cool, color=PCOL["Dry"], alpha=0.85, zorder=0)
    ax.axvspan(d_cool, t_d[-1], color=PCOL["Cooldown"], alpha=0.85, zorder=0)
    ax.axvline(d_wu, color="k", ls=":", lw=0.7, alpha=0.5, zorder=1)
    ax.axvline(d_dry, color="k", ls=":", lw=0.7, alpha=0.5, zorder=1)
    ax.axvline(d_cool, color="k", ls=":", lw=0.7, alpha=0.5, zorder=1)
    # Major grid every 5 d; minor every 1 d
    ax.xaxis.set_major_locator(mticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax.grid(which="major", axis="both", ls="-", lw=0.5, alpha=0.40, zorder=1)
    ax.grid(which="minor", axis="x", ls=":", lw=0.3, alpha=0.25, zorder=1)


# ── Figure 1: time series (5 subplots) ───────────────────────────────────────
fig1, axs = plt.subplots(5, 1, figsize=(11, 12), sharex=True)
fig1.suptitle("GWSWEX demo – loam column GWT transition (implicit solver)", fontsize=11)

_tr = mtransforms.blended_transform_factory(axs[0].transData, axs[0].transAxes)
for lbl, xd in [
    ("Warmup", d_wu / 2),
    ("Wet", (d_wu + d_dry) / 2),
    ("Dry", (d_dry + d_cool) / 2),
    ("Cooldown", (d_cool + t_d[-1]) / 2),
]:
    axs[0].text(xd, 0.97, lbl, ha="center", va="top", fontsize=8, color="#444", transform=_tr)

# Panel 0: GW head
ax = axs[0]
_shade(ax)
ax.plot(t_d, np.minimum(ts_gwh, Z_TOP), color="#882255", lw=1.8, zorder=2, label="GW head")
ax.axhline(Z_TOP, color="#8b3a07", ls="--", lw=1.0, zorder=2, label=f"Surface ({Z_TOP} m)")
ax.set_ylabel("GW head [m]")
ax.legend(fontsize=8, loc="lower right")

# Panel 1: Surface water
ax = axs[1]
_shade(ax)
ax.fill_between(t_d, ts_sw * 1e3, 0, color="#4393c3", alpha=0.6, zorder=2)
ax.plot(t_d, ts_sw * 1e3, color="#1565a0", lw=1.0, zorder=2)
ax.set_ylabel("SW depth [mm]")

# Panel 2: Column-average saturation + total UZ storage (twin y)
ax = axs[2]
_shade(ax)
sat_data = ts_sat
ax.set_ylim(np.min(sat_data) - 0.1, max(1.0, np.max(sat_data) + 0.1))
ax.plot(t_d, sat_data, color="#1a7a1a", lw=1.8, zorder=2, label=r"$\bar{\theta}$ (left)")
ax.set_ylabel(r"Column avg $\bar{\theta}$ [--]", color="#1a7a1a")
ax.tick_params(axis="y", labelcolor="#1a7a1a")
ax2_uz = ax.twinx()
ax2_uz.plot(t_d, ts_uz * 1e3, color="#8b5e00", lw=1.5, ls="--", zorder=2, label="UZ total (right)")
ax2_uz.set_ylabel("Total UZ storage [mm]", color="#8b5e00")
ax2_uz.tick_params(axis="y", labelcolor="#8b5e00")
# Combined legend
lines1, labs1 = ax.get_legend_handles_labels()
lines2, labs2 = ax2_uz.get_legend_handles_labels()
ax.legend(lines1 + lines2, labs1 + labs2, fontsize=7, loc="upper right")

# Panel 3: Cumulative fluxes -- P, PE+PT (potential), AE+AT (actual)
ax = axs[3]
_shade(ax)
ax.plot(t_d, cum_ae, color="#a50026", lw=1.5, ls="-", zorder=2, label="Cum AE (actual)")
ax.plot(t_d, cum_at, color="#d73027", lw=1.5, ls=":", zorder=2, label="Cum AT (actual)")
ax.plot(t_d, cum_ae + cum_at, color="#7f2704", lw=1.8, ls="-", zorder=2, label="Cum AET (actual)")
ax.set_ylabel("Cumulative flux [mm]")
ax.legend(fontsize=7, loc="upper left", ncol=2)

# Panel 4: Cumulative intrinsic MB error
ax = axs[4]
_shade(ax)
ax.plot(t_d, mb_cum, color="#7b3f7f", lw=1.6, zorder=2)
ax.axhline(0, color="k", lw=0.6, alpha=0.5, zorder=2)
ax.set_ylabel("Cum. intrinsic\nMB error [mm]")
ax.set_xlabel("Time [days]")

fig1.tight_layout()
fig1.savefig(OUTDIR / "timeseries.png", dpi=150, bbox_inches="tight")
print(f"Saved: {OUTDIR / 'timeseries.png'}")
plt.close(fig1)

# ── Figure 2: soil saturation profiles at key times ──────────────────────────
n_snaps = len(SNAP_ORDER)
fig2, axs2 = plt.subplots(1, n_snaps, figsize=(3.0 * n_snaps, 5.5), sharey=True)
fig2.suptitle(
    r"Soil moisture profiles  $\theta$  [m$^3$ m$^{-3}$]  at key times",
    fontsize=11,
)

for ax2, lbl in zip(axs2, SNAP_ORDER):
    theta_v = snap_theta[lbl]
    gwh = snap_gwh[lbl]
    gw_depth = max(0.0, Z_TOP - gwh)

    x_step = np.repeat(theta_v, 2)
    y_step = np.empty(2 * NL)
    y_step[0::2] = DEPTH_BNDS[:-1]
    y_step[1::2] = DEPTH_BNDS[1:]
    ax2.plot(x_step, y_step, color="#882255", lw=1.8)

    for db in DEPTH_BNDS[1:-1]:
        ax2.axhline(db, color="grey", ls=":", lw=0.5, alpha=0.55)

    # Reference lines for theta_r and theta_s
    ax2.axvline(THETA_R, color="#888", ls=":", lw=0.8, alpha=0.7)
    ax2.axvline(THETA_S, color="#888", ls=":", lw=0.8, alpha=0.7)
    ax2.axhline(gw_depth, color="#b2182b", ls="--", lw=1.2, label=f"GWT {gw_depth:.2f} m")
    ax2.legend(fontsize=6.5, loc="lower right")
    ax2.set_xlim(THETA_R - 0.01, THETA_S + 0.02)
    ax2.set_xlabel(r"$\theta$ [m$^3$ m$^{-3}$]", fontsize=9)
    ax2.set_title(lbl, fontsize=8.5)
    ax2.grid(which="major", ls=":", lw=0.4, alpha=0.45)
    ax2.xaxis.set_major_locator(mticker.MultipleLocator(0.1))
    ax2.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax2.yaxis.set_minor_locator(mticker.MultipleLocator(0.1))

axs2[0].invert_yaxis()
axs2[0].set_ylabel("Depth below surface [m]")

fig2.tight_layout()
fig2.savefig(OUTDIR / "profiles.png", dpi=150, bbox_inches="tight")
print(f"Saved: {OUTDIR / 'profiles.png'}")
plt.close(fig2)

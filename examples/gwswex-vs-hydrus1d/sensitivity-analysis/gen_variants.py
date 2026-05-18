"""Generate per-soil variants of the comparison notebooks.

Reads:
    examples/gwswex-vs-hydrus1d/comparison.ipynb              (basic, loam)
    examples/gwswex-vs-hydrus1d/comparison-intensive.ipynb    (intensive, loam)

Writes for each (setup x soil):
    comparison{-intensive}-{soil_tag}.ipynb

Output paths inside the notebooks are rewritten from
    outputs/{basic|intensive}/...
to
    outputs/{basic|intensive}/{soil_tag}/...

Existing loam notebooks are left in place; the generator writes
comparison-loam.ipynb and comparison-intensive-loam.ipynb as the
canonical loam variants. Removal of the original two is done by the
caller.

Approach: surgical text replacement on each cell's source.
"""

from __future__ import annotations

import json
import re
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EX_DIR = ROOT / "examples" / "gwswex-vs-hydrus1d"
# Source templates: the loam variants. After the basic-rename pass the basic
# template lives at `comparison-basic-loam.ipynb`; fall back to the legacy name.
_SRC_BASIC_NEW = EX_DIR / "comparison-basic-loam.ipynb"
_SRC_BASIC_OLD = EX_DIR / "comparison-loam.ipynb"
SRC_BASIC = _SRC_BASIC_NEW if _SRC_BASIC_NEW.exists() else _SRC_BASIC_OLD
SRC_INT = EX_DIR / "comparison-intensive-loam.ipynb"

# ---------------------------------------------------------------------------
# Single source of truth for soils, geometry, schedule, vegetation, forcings.
# All variant-fixed (non-OAT-tuned) experiment definitions live in JSON so
# they can be edited without touching this generator.
# ---------------------------------------------------------------------------
SOT = json.loads((EX_DIR / "experiment_definitions.json").read_text())
SOIL_DB = {name: {k: v for k, v in d.items() if not k.startswith("_")} for name, d in SOT["soils"].items()}

# Internal VARIANTS view used by main(): per soil_tag -> (LAYERS, basic_forcing, intensive_forcing)
# Built from the JSON SoT.
VARIANTS = {
    soil_tag: (
        [tuple(layer) for layer in v["layers"]],
        dict(v["basic"]),
        dict(v["intensive"]),
    )
    for soil_tag, v in SOT["variants"].items()
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _src(cell) -> str:
    s = cell["source"]
    return "".join(s) if isinstance(s, list) else s


def _set_src(cell, text: str) -> None:
    # Notebook cells store source as list of lines (with trailing newlines on
    # all but the last). Use splitlines(keepends=True) which preserves them.
    cell["source"] = text.splitlines(keepends=True)


def _build_materials_block(layers) -> str:
    """Produce a Python string assigning SOIL_TAG, LAYERS, MATERIALS, LAYER_MAT_IDS, THETA_R, THETA_S."""
    return (
        "# ── Soil property catalogue (Carsel & Parrish 1988 defaults) ─────────────────\n"
        "SOIL_DB = {\n"
        '    "sand": dict(theta_r=0.045, theta_s=0.430, alpha=14.5, n=2.68, K_sat=7.128,  lam=0.5),\n'
        '    "loam": dict(theta_r=0.078, theta_s=0.430, alpha=3.6,  n=1.56, K_sat=0.2496, lam=0.5),\n'
        '    "clay": dict(theta_r=0.068, theta_s=0.380, alpha=0.8,  n=1.09, K_sat=0.0480, lam=0.5),\n'
        "}\n"
        "\n"
        f"SOIL_TAG = {layers!r}\n"  # placeholder; replaced below
    )


def _layers_python_repr(layers) -> str:
    """Render LAYERS list using fraction-form for readable diffs."""
    parts = []
    for name, frac in layers:
        # Try to render thirds nicely
        if abs(frac - 1.0) < 1e-12:
            fr = "1.0"
        elif abs(frac - 1 / 3) < 1e-9:
            fr = "1/3"
        elif abs(frac - 2 / 3) < 1e-9:
            fr = "2/3"
        elif abs(frac - 0.5) < 1e-12:
            fr = "0.5"
        else:
            fr = repr(frac)
        parts.append(f'("{name}", {fr})')
    return "[" + ", ".join(parts) + "]"


def _soil_db_block_lines() -> list[str]:
    """Render SOIL_DB literally from the JSON SoT (single source of truth).

    Each soil row is dropped with a brief inline annotation if relaxed away
    from Carsel & Parrish defaults (e.g. clay).
    """
    lines = [
        "# ── Soil property catalogue (loaded from experiment_definitions.json) ────────",
        "# Single source of truth: examples/gwswex-vs-hydrus1d/experiment_definitions.json",
        "SOIL_DB = {",
    ]
    for name, props in SOIL_DB.items():
        src = SOT["soils"][name].get("_source", "")
        comment = f"  # {src}" if src else ""
        lines.append(
            f'    "{name}": dict(theta_r={props["theta_r"]:.3f}, theta_s={props["theta_s"]:.3f}, '
            f'alpha={props["alpha"]:g}, n={props["n"]:g}, K_sat={props["K_sat"]:g}, '
            f'lam={props["lam"]:g}),{comment}'
        )
    lines.append("}")
    return lines


def make_soil_block(soil_tag: str, layers) -> str:
    return (
        "\n".join(_soil_db_block_lines()) + "\n\n"
        f'SOIL_TAG = "{soil_tag}"\n'
        f"LAYERS = {_layers_python_repr(layers)}    # (soil_name, fraction_of_column_top_to_bottom)\n"
        "\n"
        "# Build materials list and per-layer assignment (top-to-bottom).\n"
        "MATERIALS = []\n"
        "_seen: dict[str, int] = {}\n"
        "for _name, _ in LAYERS:\n"
        "    if _name in _seen:\n"
        "        continue\n"
        "    _seen[_name] = len(MATERIALS) + 1\n"
        "    _s = SOIL_DB[_name]\n"
        "    MATERIALS.append(dict(\n"
        '        id=_seen[_name], name=_name, K_sat=_s["K_sat"], lam=_s["lam"],\n'
        '        vanG=dict(alpha=_s["alpha"], n=_s["n"], theta_r=_s["theta_r"], theta_s=_s["theta_s"]),\n'
        "    ))\n"
        "\n"
        "LAYER_MAT_IDS: list[int] = []\n"
        "for _name, _frac in LAYERS:\n"
        "    LAYER_MAT_IDS += [_seen[_name]] * int(round(NL * _frac))\n"
        "# Pad/trim to exactly NL layers (handles fractional rounding).\n"
        "LAYER_MAT_IDS = (LAYER_MAT_IDS + [LAYER_MAT_IDS[-1]] * NL)[:NL]\n"
        "\n"
        "# Bottom-material θr / θs: GW zone resides here. Used for HYDRUS plotting\n"
        "# scalar references and the residual saturated-zone correction in the MB closure.\n"
        "_bot_soil = SOIL_DB[LAYERS[-1][0]]\n"
        'THETA_R, THETA_S = _bot_soil["theta_r"], _bot_soil["theta_s"]\n'
        "\n"
    )


def make_hydrus_material_block() -> str:
    """Replacement for the single-material HYDRUS material+profile setup block."""
    return (
        "mat_df = ml.get_empty_material_df(n=len(MATERIALS))\n"
        "for _mat in MATERIALS:\n"
        '    _s = SOIL_DB[_mat["name"]]\n'
        '    mat_df.loc[_mat["id"]] = [_s["theta_r"], _s["theta_s"], _s["alpha"]/100.0, _s["n"], _s["K_sat"]*100.0, _s["lam"]]\n'
        "ml.add_material(mat_df)\n"
    )


_DEMO_COMMENT_PATTERNS = [
    # Strip lines or in-line tails referring to demo configuration.
    (re.compile(r"\s*\(demo configuration\)"), ""),
    (re.compile(r"\s*\(matches demo plotting\)"), ""),
    (re.compile(r"\s*\(matches demo\)"), ""),
    (re.compile(r"\s*\(no growth — matches demo\)"), ""),
    (re.compile(r"\s*— match demo$", re.MULTILINE), ""),
    (re.compile(r" identical to the demos\)"), ")"),
    (re.compile(r" identical to the demos:"), ":"),
    (re.compile(r" — same as the demo[^\n]*"), ""),
    (
        re.compile(
            r"^.*Geometry, vegetation, IC, and forcing schedule match `demo-explicit` /\n.*`demo-implicit` exactly\.[^\n]*\n",
            re.MULTILINE,
        ),
        "",
    ),
    (re.compile(r" identical to the demos\.", re.IGNORECASE), "."),
    (re.compile(r" identical to the demos", re.IGNORECASE), ""),
    (re.compile(r" used in demo-explicit\.py / demo-implicit\.py"), ""),
    (
        re.compile(r"# ── Vegetation: pasture, fixed roots at 0\.6 m \(demo configuration\) ───────────"),
        "# ── Vegetation: pasture, fixed roots at 0.6 m ────────────────────────────────",
    ),
    (
        re.compile(r"# ── GWSWEX physics parameters \(matches demo\) ─────────────────────────────────"),
        "# ── GWSWEX physics parameters ────────────────────────────────────────────────",
    ),
    (
        re.compile(r"# Constant root depth at 60 cm \(no growth — matches demo\)"),
        "# Constant root depth at 60 cm (no growth)",
    ),
    (re.compile(r"# Phase shading \(matches demo plotting\)"), "# Phase shading"),
    (re.compile(r"# Snapshot points \(hour index → label\) — match demo"), "# Snapshot points (hour index → label)"),
]


def _strip_demo_comments(src: str) -> str:
    out = src
    for pat, rep in _DEMO_COMMENT_PATTERNS:
        out = pat.sub(rep, out)
    return out


_CLIP_BLOCK = (
    "h_gw_d = np.interp(t_d, h_t, h_gw_all)\n"
    "# Clip HYDRUS-1D WT to model top (surface): plotted depth ≥ 0 cm.\n"
    "h_gw_d = np.maximum(h_gw_d, 0.0)\n"
    "h_gw_all = np.maximum(h_gw_all, 0.0)"
)

# Idempotency: collapse any prior `np.maximum(h_gw_d|h_gw_all, 0.0)` lines and
# re-emit the canonical block from the seed `np.interp` line. This lets the
# generator be re-run on its own output without accumulating duplicates.
_CLIP_DUP_RE = re.compile(
    r"h_gw_d = np\.interp\(t_d, h_t, h_gw_all\)\s*\n"
    r"(?:[ \t]*#[^\n]*\n|"
    r"[ \t]*h_gw_d = np\.maximum\(h_gw_d, 0\.0\)\s*\n|"
    r"[ \t]*h_gw_all = np\.maximum\(h_gw_all, 0\.0\)\s*\n)+"
)


_FORCING_TWIN_DEF_RE = re.compile(
    r"def _forcing_twin\(ax\):.*?(?=\ndef |\nclass |\n# |\Z)",
    re.DOTALL,
)
_FORCING_TWIN_CALL_RE = re.compile(r"^\s*_forcing_twin\([^)]*\)\s*\n", re.MULTILINE)


def _apply_post_processing(cells, soil_tag: str, setup: str) -> None:
    """Apply across-cell mechanical fixes: clip HYDRUS WT, retitle plots, strip demo refs.

    setup: 'basic' or 'intensive'.
    """
    for c in cells:
        if c.get("cell_type") != "code":
            continue
        s = _src(c)
        # Strip demo-referencing comments
        s = _strip_demo_comments(s)
        # Strip the precip / PE+PT twin-axis overlay everywhere
        # (helper definition + every call site). Keeps phase shading intact.
        s = _FORCING_TWIN_DEF_RE.sub("", s)
        s = _FORCING_TWIN_CALL_RE.sub("", s)
        # Drop the misleading "(negative ⇒ ponded)" axis annotation: the GWSWEX
        # WT trace is now clipped at the surface and ponding is shown as its
        # own dedicated SW subplot.
        s = s.replace(
            'ax_wt.set_ylabel("WT depth [cm]   (negative ⇒ ponded)")',
            'ax_wt.set_ylabel("WT depth [cm]")',
        )
        # Clip HYDRUS WT to surface (idempotent: collapse any prior block, re-emit).
        s = _CLIP_DUP_RE.sub(_CLIP_BLOCK + "\n", s)
        # WT-difference subplot: clip BOTH HYDRUS and GWSWEX to surface so the
        # difference reflects only physically resolvable WT depths (above-surface
        # values are conventionally pinned to 0 = model top in both models).
        s = s.replace(
            'resid = h_gw_d - res["gw_depth_cm"]',
            'resid = h_gw_d - np.maximum(res["gw_depth_cm"], 0.0)',
        )
        # WT plot lines: clip GWSWEX to surface in the plot too (HYDRUS already clipped).
        s = s.replace(
            'ax_wt.plot(t_d, res_i["gw_depth_cm"], **STYLE["GWSWEX implicit"])',
            'ax_wt.plot(t_d, np.maximum(res_i["gw_depth_cm"], 0.0), **STYLE["GWSWEX implicit"])',
        )
        s = s.replace(
            'ax_wt.plot(t_d, res_e["gw_depth_cm"], **STYLE["GWSWEX explicit"])',
            'ax_wt.plot(t_d, np.maximum(res_e["gw_depth_cm"], 0.0), **STYLE["GWSWEX explicit"])',
        )
        # Per-snapshot GWSWEX axhline (depth ≥ 0).
        s = s.replace(
            'ax.axhline(res["gw_depth_cm"][d], color=STYLE[key]["color"], ls="--", lw=0.9, alpha=0.5)',
            'ax.axhline(max(res["gw_depth_cm"][d], 0.0), color=STYLE[key]["color"], ls="--", lw=0.9, alpha=0.5)',
        )
        # Metric calls: pass clipped GWSWEX to the WT-depth metrics (HYDRUS already clipped).
        for _fn in ("_rmse", "_mae", "_nse", "_bias"):
            s = s.replace(
                f"{_fn}(res['gw_depth_cm'], h_gw_d)",
                f"{_fn}(np.maximum(res['gw_depth_cm'], 0.0), h_gw_d)",
            )
        # Phase-mask metrics: clip the GWSWEX side fed to (sim, obs).
        s = s.replace(
            'sim, obs = res["gw_depth_cm"][mask], h_gw_d[mask]',
            'sim, obs = np.maximum(res["gw_depth_cm"], 0.0)[mask], h_gw_d[mask]',
        )
        # Plot title: use SOIL_TAG instead of hard-coded 'loam'
        s = s.replace(
            'fig.suptitle("Water-table depth: 65-day loam column", fontsize=11)',
            'fig.suptitle(f"Water-table depth: 65-day {SOIL_TAG} column", fontsize=11)',
        )
        s = s.replace(
            'fig.suptitle(f"Water-table depth: {T_TOTAL}-day loam column (32 d / 4-phase / cooldown)", fontsize=11)',
            'fig.suptitle(f"Water-table depth: {T_TOTAL}-day {SOIL_TAG} column (32 d / 4-phase / cooldown)", fontsize=11)',
        )
        # Typo fix carried over from loam template.
        s = s.replace('ax_res.set_ylabel("Diffeerence")', 'ax_res.set_ylabel("Difference [cm]")')

        # Zone-averaged theta plot: drop sharey=True so each panel autoscales,
        # then add a dynamic ylim per subplot (combined HYDRUS+GWSWEX trace
        # range padded by 0.05 m3/m3) to avoid clipping the deeper zones.
        s = s.replace(
            "fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True, sharey=True)",
            "fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True)",
        )
        # Inject ylim setter inside the per-zone loop, just after `_phase_lines(ax)`.
        _OLD_AX = (
            "    ax.plot(t_d, res_e[\"theta\"][:, gm].mean(1), **STYLE[\"GWSWEX explicit\"])\n"
            "    ax.set_title(label, fontsize=9)\n"
            "    _phase_lines(ax)\n"
        )
        _NEW_AX = (
            "    ax.plot(t_d, res_e[\"theta\"][:, gm].mean(1), **STYLE[\"GWSWEX explicit\"])\n"
            "    # Dynamic per-subplot ylim: combined trace range +/- 0.05 m3/m3 padding.\n"
            "    _vals = np.concatenate([h_zone, res_i[\"theta\"][:, gm].mean(1), res_e[\"theta\"][:, gm].mean(1)])\n"
            "    _ymin, _ymax = float(np.nanmin(_vals)), float(np.nanmax(_vals))\n"
            "    ax.set_ylim(max(0.0, _ymin - 0.05), min(1.0, _ymax + 0.05))\n"
            "    ax.set_title(label, fontsize=9)\n"
            "    _phase_lines(ax)\n"
        )
        if _OLD_AX in s:
            s = s.replace(_OLD_AX, _NEW_AX)
        _set_src(c, s)

    # HYDRUS solver tolerance: relax for layered intensive profiles where
    # sharp Ks contrasts at the interface stress Picard convergence.
    if setup == "intensive":
        for c in cells:
            if c.get("cell_type") != "code":
                continue
            s = _src(c)
            s = s.replace(
                "                 maxit=50, tolh=0.1, tolth=1e-4, ha=1e-6, hb=1e4)",
                "                 maxit=200, tolh=1.0, tolth=1e-3, ha=1e-6, hb=1e4)",
            )
            _set_src(c, s)

        # ---- Intensive: convert 2-panel WT|Diff figure into 3-panel WT|SW|Diff ----
        # HYDRUS-1D in the intensive setup uses top_bc=2 (atmospheric BC with
        # surface layer), so HYDRUS DOES represent ponding via h_sw_d.
        _OLD_WT_INT = (
            'fig, (ax_wt, ax_res) = plt.subplots(2, 1, figsize=(10, 5.8), '
            'gridspec_kw={"height_ratios": [3, 1]}, sharex=True)\n'
        )
        _NEW_WT_INT = (
            "fig, (ax_wt, ax_sw, ax_res) = plt.subplots(\n"
            "    3, 1, figsize=(10, 7.4),\n"
            '    gridspec_kw={"height_ratios": [3, 1, 1]}, sharex=True,\n'
            ")\n"
        )
        _OLD_RES_INT = (
            "for key, res in [(\"GWSWEX implicit\", res_i), (\"GWSWEX explicit\", res_e)]:\n"
            "    # Positive = GWSWEX WT shallower (higher head) than HYDRUS-1D\n"
            "    resid = h_gw_d - np.maximum(res[\"gw_depth_cm\"], 0.0)\n"
        )
        _NEW_SW_INT = (
            "# ── Ponded surface-water depth ───────────────────────────────────────────\n"
            "# Both models accumulate ponding (HYDRUS: top_bc=2 surface layer; GWSWEX:\n"
            "# SW reservoir). WT and SW are separate state variables: GW can sit deep\n"
            "# below the surface while a thin pond persists, so the two panels carry\n"
            "# independent information.\n"
            "ax_sw.plot(t_d, h_sw_d, **STYLE[\"HYDRUS-1D\"])\n"
            "for key, res in [(\"GWSWEX implicit\", res_i), (\"GWSWEX explicit\", res_e)]:\n"
            "    ax_sw.plot(t_d, res[\"sw_cm\"], lw=1.4, color=STYLE[key][\"color\"],\n"
            "               ls=STYLE[key][\"ls\"], label=key)\n"
            "_sw_top = max(float(h_sw_d.max()), float(res_i[\"sw_cm\"].max()),\n"
            "              float(res_e[\"sw_cm\"].max()), 1e-3)\n"
            "ax_sw.fill_between(t_d, np.maximum.reduce([h_sw_d, res_i[\"sw_cm\"], res_e[\"sw_cm\"]]),\n"
            "                   0, alpha=0.10, color=\"steelblue\")\n"
            "ax_sw.set_ylabel(\"Ponded SW [cm]\")\n"
            "ax_sw.axhline(0, color=\"k\", lw=0.5, alpha=0.5)\n"
            "ax_sw.set_ylim(-0.05 * _sw_top, 1.1 * _sw_top)\n"
            "ax_sw.legend(loc=\"upper right\", fontsize=8)\n"
            "_shade_phases(ax_sw, label=False)\n"
            "\n"
            "for key, res in [(\"GWSWEX implicit\", res_i), (\"GWSWEX explicit\", res_e)]:\n"
            "    # Positive = GWSWEX WT shallower (higher head) than HYDRUS-1D\n"
            "    resid = h_gw_d - np.maximum(res[\"gw_depth_cm\"], 0.0)\n"
        )
        for c in cells:
            if c.get("cell_type") != "code":
                continue
            s = _src(c)
            if _OLD_WT_INT in s and _OLD_RES_INT in s:
                s = s.replace(_OLD_WT_INT, _NEW_WT_INT)
                s = s.replace(_OLD_RES_INT, _NEW_SW_INT)
                _set_src(c, s)

    # Basic plot harmonisation: replace the 2-line `_phase_lines` definition
    # with a 3-phase shaded `PHASES` + `_shade_phases` block (matching the
    # intensive notebook), and migrate all `_phase_lines(...)` call sites.
    if setup == "basic":
        _OLD_DEF = (
            "def _phase_lines(ax):\n"
            "    ylo, yhi = ax.get_ylim()\n"
            "    for t, label in [(T_P2, \"Wet\"), (T_P3, \"Dry\")]:\n"
            "        ax.axvline(t, color=\"k\", ls=\":\", lw=0.8, alpha=0.6)\n"
            "        ax.text(t + 0.4, yhi - 0.03 * (yhi - ylo), label, fontsize=8, va=\"top\", alpha=0.65)"
        )
        _NEW_DEF = (
            "# Phase shading (Warmup / Wet / Dry) - harmonised with intensive notebook layout.\n"
            "PHASES = [\n"
            "    (T_P1, T_P2, \"Warmup\", \"#ddeef7\"),\n"
            "    (T_P2, T_P3, \"Wet\",    \"#ceecc8\"),\n"
            "    (T_P3, t_d[-1], \"Dry\", \"#faecd1\"),\n"
            "]\n"
            "\n"
            "def _shade_phases(ax, label=True):\n"
            "    ylo, yhi = ax.get_ylim()\n"
            "    for x0, x1, lbl, col in PHASES:\n"
            "        ax.axvspan(x0, x1, color=col, alpha=0.55, zorder=0)\n"
            "        if label:\n"
            "            ax.text(0.5 * (x0 + x1), yhi - 0.04 * (yhi - ylo),\n"
            "                    lbl, fontsize=8, ha=\"center\", va=\"top\", alpha=0.75)\n"
            "    ax.set_ylim(ylo, yhi)\n"
            "\n"
            "def _phase_lines(ax, label=True):\n"
            "    \"\"\"Backward-compat shim - delegates to _shade_phases.\"\"\"\n"
            "    _shade_phases(ax, label=label)"
        )
        for c in cells:
            if c.get("cell_type") != "code":
                continue
            s = _src(c)
            if _OLD_DEF in s:
                s = s.replace(_OLD_DEF, _NEW_DEF)
                _set_src(c, s)

        # ---- Basic-only: capture ponded SW depth and add SW subplot to WT figure ----
        # 1) Extend run_and_collect to record state["SW"][0] each step and
        #    expose `sw_cm` in the result dict. Idempotent: each substitution
        #    only fires when the new form is not already present.
        for c in cells:
            if c.get("cell_type") != "code":
                continue
            s = _src(c)
            if "sw_all" not in s:
                s = s.replace(
                    "    gw_all, theta_all, cum_e_all, cum_t_all, mbi_all = [], [], [], [], []",
                    "    gw_all, sw_all, theta_all, cum_e_all, cum_t_all, mbi_all = [], [], [], [], [], []",
                )
                s = s.replace(
                    '        gw_all.append(state["GWH"][0])\n'
                    '        theta_all.append(state["theta"][:, 0].copy())\n',
                    '        gw_all.append(state["GWH"][0])\n'
                    '        sw_all.append(float(state["SW"][0]))\n'
                    '        theta_all.append(state["theta"][:, 0].copy())\n',
                )
            if "sw_cm=sw_arr" not in s:
                s = s.replace(
                    "    return dict(\n" "        gw_depth_cm=(Z_TOP - gw_arr[idx]) * 100,\n",
                    "    sw_arr = np.array(sw_all)\n"
                    "    return dict(\n"
                    "        gw_depth_cm=(Z_TOP - gw_arr[idx]) * 100,\n"
                    "        sw_cm=sw_arr[idx] * 100.0,\n",
                )
            _set_src(c, s)

        # 2) Three-panel WT figure: WT (3) | SW (1) | Difference (1).
        _OLD_WT = (
            'fig, (ax_wt, ax_res) = plt.subplots(2, 1, figsize=(9, 5.5), '
            'gridspec_kw={"height_ratios": [3, 1]}, sharex=True)\n'
        )
        _NEW_WT = (
            "fig, (ax_wt, ax_sw, ax_res) = plt.subplots(\n"
            "    3, 1, figsize=(9, 7.0),\n"
            '    gridspec_kw={"height_ratios": [3, 1, 1]}, sharex=True,\n'
            ")\n"
        )
        _OLD_RES_BLOCK = (
            "for key, res in [(\"GWSWEX implicit\", res_i), (\"GWSWEX explicit\", res_e)]:\n"
            "    # positive = GWSWEX WT is shallower (higher head) than HYDRUS-1D\n"
            "    # negative = GWSWEX WT is deeper (lower head) than HYDRUS-1D\n"
            "    resid = h_gw_d - np.maximum(res[\"gw_depth_cm\"], 0.0)\n"
        )
        _NEW_SW_BLOCK = (
            "# ── Ponded surface-water depth (GWSWEX SW reservoir) ─────────────────────\n"
            "# HYDRUS-1D in the basic setup uses top_bc=3 (atmospheric BC with surface\n"
            "# runoff): excess infiltration is shed, so HYDRUS does not represent ponding.\n"
            "# Only GWSWEX traces are shown.\n"
            "for key, res in [(\"GWSWEX implicit\", res_i), (\"GWSWEX explicit\", res_e)]:\n"
            "    ax_sw.plot(t_d, res[\"sw_cm\"], lw=1.4, color=STYLE[key][\"color\"],\n"
            "               ls=STYLE[key][\"ls\"], label=key)\n"
            "ax_sw.fill_between(t_d, np.maximum(res_i[\"sw_cm\"], res_e[\"sw_cm\"]), 0,\n"
            "                   alpha=0.10, color=\"steelblue\")\n"
            "ax_sw.set_ylabel(\"Ponded SW [cm]\")\n"
            "ax_sw.axhline(0, color=\"k\", lw=0.5, alpha=0.5)\n"
            "_sw_max = max(float(res_i[\"sw_cm\"].max()), float(res_e[\"sw_cm\"].max()), 1e-3)\n"
            "ax_sw.set_ylim(-0.05 * _sw_max, 1.1 * _sw_max)\n"
            "ax_sw.legend(loc=\"upper right\", fontsize=8)\n"
            "_phase_lines(ax_sw, label=False)\n"
            "\n"
            "for key, res in [(\"GWSWEX implicit\", res_i), (\"GWSWEX explicit\", res_e)]:\n"
            "    # positive = GWSWEX WT is shallower (higher head) than HYDRUS-1D\n"
            "    # negative = GWSWEX WT is deeper (lower head) than HYDRUS-1D\n"
            "    resid = h_gw_d - np.maximum(res[\"gw_depth_cm\"], 0.0)\n"
        )
        for c in cells:
            if c.get("cell_type") != "code":
                continue
            s = _src(c)
            if _OLD_WT in s and _OLD_RES_BLOCK in s:
                s = s.replace(_OLD_WT, _NEW_WT)
                s = s.replace(_OLD_RES_BLOCK, _NEW_SW_BLOCK)
                _set_src(c, s)


def make_hydrus_profile_per_node_mat_block(z_col_cm_expr: str) -> str:
    """Snippet to assign per-node Mat ID after the profile is created."""
    return (
        "# Per-node material assignment based on depth (top-to-bottom layered profile).\n"
        "_z_node_depth_cm = -profile[\"x\"].to_numpy(float)        # 0..Z_COL_CM, 0=surface\n"
        f"_dz_layer_cm = ({z_col_cm_expr}) / NL\n"
        "_node_layer_idx = np.minimum(np.floor(_z_node_depth_cm / _dz_layer_cm).astype(int), NL - 1)\n"
        "profile[\"Mat\"] = np.array([LAYER_MAT_IDS[i] for i in _node_layer_idx], dtype=int)\n"
    )


# ---------------------------------------------------------------------------
# Cell-2 builders — load fixed quantities from experiment_definitions.json
# (single source of truth). Only OAT-tuned MODEL_PARAMS / ET_STRESS remain
# literal so apply_oat_optima can rewrite them in-place via regex.
# ---------------------------------------------------------------------------

_CELL2_SOT_PRELUDE = (
    "# ── Single source of truth: experiment_definitions.json ──────────────────────\n"
    "# Soils, geometry, schedule, vegetation, layering and per-variant atmospheric\n"
    "# forcings are loaded from the SoT JSON sitting next to this notebook. Only\n"
    "# tuned MODEL_PARAMS / ET_STRESS / set_solver kwargs (sensitivity-analysis\n"
    "# outputs) remain literal in this notebook.\n"
    "import json\n"
    "_SOT = json.loads((HERE / \"experiment_definitions.json\").read_text())\n"
)

_CELL2_SOIL_DERIV = (
    "# Soils + variant from SoT\n"
    "SOIL_DB = {n: {k: v for k, v in d.items() if not k.startswith(\"_\")}\n"
    "           for n, d in _SOT[\"soils\"].items()}\n"
    "_var = _SOT[\"variants\"][SOIL_TAG]\n"
    "LAYERS = [tuple(_l) for _l in _var[\"layers\"]]    # (soil_name, fraction_top_to_bottom)\n"
    "\n"
    "# Build materials list and per-layer assignment (top-to-bottom).\n"
    "MATERIALS = []\n"
    "_seen: dict[str, int] = {}\n"
    "for _name, _ in LAYERS:\n"
    "    if _name in _seen:\n"
    "        continue\n"
    "    _seen[_name] = len(MATERIALS) + 1\n"
    "    _s = SOIL_DB[_name]\n"
    "    MATERIALS.append(dict(\n"
    "        id=_seen[_name], name=_name, K_sat=_s[\"K_sat\"], lam=_s[\"lam\"],\n"
    "        vanG=dict(alpha=_s[\"alpha\"], n=_s[\"n\"], theta_r=_s[\"theta_r\"], theta_s=_s[\"theta_s\"]),\n"
    "    ))\n"
    "\n"
    "LAYER_MAT_IDS: list[int] = []\n"
    "for _name, _frac in LAYERS:\n"
    "    LAYER_MAT_IDS += [_seen[_name]] * int(round(NL * _frac))\n"
    "# Pad/trim to exactly NL layers (handles fractional rounding).\n"
    "LAYER_MAT_IDS = (LAYER_MAT_IDS + [LAYER_MAT_IDS[-1]] * NL)[:NL]\n"
    "\n"
    "# Bottom-material θr / θs: GW zone resides here.\n"
    "_bot_soil = SOIL_DB[LAYERS[-1][0]]\n"
    "THETA_R, THETA_S = _bot_soil[\"theta_r\"], _bot_soil[\"theta_s\"]\n"
    "\n"
)


def _make_cell2_basic(soil_tag: str) -> str:
    return (
        _CELL2_SOT_PRELUDE
        + f'SETUP, SOIL_TAG = "basic", "{soil_tag}"\n'
        + "\n"
        + "# ── Geometry [m] ─────────────────────────────────────────────────────────────\n"
        + "_g = _SOT[\"geometry\"][SETUP]\n"
        + "Z_TOP, Z_BOT, DZ = _g[\"Z_TOP\"], _g[\"Z_BOT\"], _g[\"DZ\"]\n"
        + "NL, Z_WT = _g[\"NL\"], _g[\"Z_WT\"]\n"
        + "NE = 1\n"
        + "BNDS = np.linspace(Z_TOP, Z_BOT, NL + 1)  # top-down, shape (NL+1,)\n"
        + "\n"
        + _CELL2_SOIL_DERIV
        + "# ── Phase structure and daily atmospheric forcing [cm d^-1] ──────────────────\n"
        + "_sch = _SOT[\"schedule\"][SETUP]\n"
        + "T_TOTAL = _sch[\"T_TOTAL_d\"]\n"
        + "T_P1, T_P2, T_P3 = _sch[\"T_P1_d\"], _sch[\"T_P2_d\"], _sch[\"T_P3_d\"]\n"
        + "\n"
        + "_f = _var[SETUP]\n"
        + "P_WET, PE_WET, PT_WET = _f[\"P_WET\"], _f[\"PE_WET\"], _f[\"PT_WET\"]   # wet-phase rates [cm d^-1]\n"
        + "PE_DRY, PT_DRY        = _f[\"PE_DRY\"], _f[\"PT_DRY\"]                # dry-phase rates [cm d^-1]\n"
        + "\n"
        + "t_d = np.arange(1, T_TOTAL + 1, dtype=float)\n"
        + "prec_d = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, P_WET, 0.0))\n"
        + "pet_d  = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, PE_WET, PE_DRY))\n"
        + "ptt_d  = np.where(t_d <= T_P2, 0.0, np.where(t_d <= T_P3, PT_WET, PT_DRY))\n"
        + "\n"
        + "# ── Model physics parameters (Green-Ampt infiltration + connectivity) ─────────\n"
        + "MODEL_PARAMS = dict(psi_f=0.01, F_min=1e-7, ICratio_min=0.20)\n"
        + "\n"
        + "# ── ET stress thresholds (Laio 2001) — attached to vegetation type ────────────\n"
        + "ET_STRESS = dict(s_star=0.4, s_w=0.1, s_h=0.05, s_e=0.3)\n"
        + "\n"
        + "# ── Root growth (linear over T_TOTAL d) — from SoT ───────────────────────────\n"
        + "_veg = _SOT[\"vegetation\"][SETUP]\n"
        + "ROOT_D0, ROOT_D1 = _veg[\"root_d_initial\"], _veg[\"root_d_final\"]  # [m]\n"
        + "T0 = datetime(2024, 1, 1)\n"
        + "\n"
        + "# ── GWSWEX layer mid-depths (top-to-bottom, cm below surface) ─────────────────\n"
        + "LAYER_DEPTH = np.linspace(DZ / 2, (Z_TOP - Z_BOT) - DZ / 2, NL) * 100.0\n"
    )


def _make_cell2_intensive(soil_tag: str) -> str:
    return (
        _CELL2_SOT_PRELUDE
        + f'SETUP, SOIL_TAG = "intensive", "{soil_tag}"\n'
        + "\n"
        + "# ── Geometry [m] ─────────────────────────────────────────────────────────────\n"
        + "_g = _SOT[\"geometry\"][SETUP]\n"
        + "Z_TOP, Z_BOT, DZ = _g[\"Z_TOP\"], _g[\"Z_BOT\"], _g[\"DZ\"]\n"
        + "NL, Z_WT = _g[\"NL\"], _g[\"Z_WT\"]\n"
        + "NE = 1\n"
        + "BNDS = np.linspace(Z_TOP, Z_BOT, NL + 1)          # top-down [m]\n"
        + "LAYER_DEPTH = np.arange(0.5, NL, 1.0)             # midpoints below surface [cm]\n"
        + "\n"
        + _CELL2_SOIL_DERIV
        + "# Convenience scalars (bottom material) — used by HYDRUS-1D head-IC and hourly K conversion.\n"
        + "K_SAT_MD = SOIL_DB[LAYERS[-1][0]][\"K_sat\"]       # bottom-material K_sat [m d^-1]\n"
        + "K_SAT_MH = K_SAT_MD / 24                          # [m h^-1]\n"
        + "LAM = SOIL_DB[LAYERS[-1][0]][\"lam\"]\n"
        + "\n"
        + "# ── Phase structure (warmup → wet → dry → cooldown) ──────────────────────────\n"
        + "_sch = _SOT[\"schedule\"][SETUP]\n"
        + "T_WU, T_WET, T_DRY, T_COOL = _sch[\"T_WU_h\"], _sch[\"T_WET_h\"], _sch[\"T_DRY_h\"], _sch[\"T_COOL_h\"]    # [h]\n"
        + "N = T_WU + T_WET + T_DRY + T_COOL                 # total hours\n"
        + "T_TOTAL = N // 24                                 # days\n"
        + "\n"
        + "t_d = (np.arange(1, N + 1)) / 24.0                # hourly axis in days\n"
        + "T_P1, T_P2, T_P3, T_P4 = 0.0, T_WU / 24.0, (T_WU + T_WET) / 24.0, (T_WU + T_WET + T_DRY) / 24.0\n"
        + "\n"
        + "# ── Atmospheric forcing ──────────────────────────────────────────────────────\n"
        + "# GWSWEX (T='h', L='m'): hourly arrays in m h^-1\n"
        + "_f = _var[SETUP]\n"
        + "P_W, E_W, T_W = _f[\"P_W\"], _f[\"E_W\"], _f[\"T_W\"]             # wet\n"
        + "E_D, T_D       = _f[\"E_D\"], _f[\"T_D\"]                       # dry\n"
        + "\n"
        + "_h_idx = np.arange(N)\n"
        + "prec_h = np.where(_h_idx < T_WU, 0.0,\n"
        + "         np.where(_h_idx < T_WU + T_WET, P_W, 0.0))\n"
        + "pet_h  = np.where(_h_idx < T_WU, 0.0,\n"
        + "         np.where(_h_idx < T_WU + T_WET, E_W,\n"
        + "         np.where(_h_idx < T_WU + T_WET + T_DRY, E_D, 0.0)))\n"
        + "ptt_h  = np.where(_h_idx < T_WU, 0.0,\n"
        + "         np.where(_h_idx < T_WU + T_WET, T_W,\n"
        + "         np.where(_h_idx < T_WU + T_WET + T_DRY, T_D, 0.0)))\n"
        + "\n"
        + "# Hourly forcing in cm h^-1 (for plotting overlay)\n"
        + "M_PER_H_TO_CM_PER_H = 100.0\n"
        + "prec_d_cmh = prec_h * M_PER_H_TO_CM_PER_H\n"
        + "pet_d_cmh  = pet_h  * M_PER_H_TO_CM_PER_H\n"
        + "ptt_d_cmh  = ptt_h  * M_PER_H_TO_CM_PER_H\n"
        + "\n"
        + "# Daily forcing arrays in cm d^-1 for HYDRUS-1D atmospheric BC.\n"
        + "M_PER_H_TO_CM_PER_D = 100.0 * 24.0\n"
        + "day_idx = np.arange(1, T_TOTAL + 1, dtype=float)\n"
        + "P_W_CMD = P_W * M_PER_H_TO_CM_PER_D\n"
        + "E_W_CMD = E_W * M_PER_H_TO_CM_PER_D\n"
        + "T_W_CMD = T_W * M_PER_H_TO_CM_PER_D\n"
        + "E_D_CMD = E_D * M_PER_H_TO_CM_PER_D\n"
        + "T_D_CMD = T_D * M_PER_H_TO_CM_PER_D\n"
        + "prec_d = np.where((day_idx >  3) & (day_idx <= 13), P_W_CMD, 0.0)\n"
        + "pet_d  = np.where((day_idx >  3) & (day_idx <= 13), E_W_CMD,\n"
        + "         np.where((day_idx > 13) & (day_idx <= 25), E_D_CMD, 0.0))\n"
        + "ptt_d  = np.where((day_idx >  3) & (day_idx <= 13), T_W_CMD,\n"
        + "         np.where((day_idx > 13) & (day_idx <= 25), T_D_CMD, 0.0))\n"
        + "\n"
        + "# ── Vegetation (from SoT) ────────────────────────────────────────────────────\n"
        + "_veg = _SOT[\"vegetation\"][SETUP]\n"
        + "ROOT_D0, ROOT_D1 = _veg[\"root_d_initial\"], _veg[\"root_d_final\"]   # [m]\n"
        + "ET_STRESS = dict(s_star=0.5, s_w=0.1, s_h=0.05, s_e=0.5)\n"
        + "\n"
        + "# ── GWSWEX physics parameters (matches demo) ─────────────────────────────────\n"
        + "MODEL_PARAMS = dict(psi_f=0.09, F_min=1e-6, ICratio_min=0.05)\n"
    )


# ---------------------------------------------------------------------------
# Basic notebook transformer
# ---------------------------------------------------------------------------


def transform_basic(orig: dict, soil_tag: str, layers, forcing: dict) -> dict:
    nb = deepcopy(orig)
    cells = nb["cells"]

    # Cell 0: markdown title — append soil tag
    md = _src(cells[0])
    md = re.sub(
        r"^# GWSWEX vs HYDRUS-1D: Solver Comparison.*$",
        f"# GWSWEX vs HYDRUS-1D: Solver Comparison — {soil_tag}",
        md,
        count=1,
        flags=re.MULTILINE,
    )
    md = re.sub(
        r"\| Column \| .*?\|",
        f"| Column | 3 m, 150 × 0.02 m layers, **{soil_tag}** |",
        md,
        count=1,
    )
    md = re.sub(
        r"Outputs \(HYDRUS-1D workspace\) are written under `outputs/basic/(?:loam/)?`\.",
        f"Outputs (HYDRUS-1D workspace) are written under `outputs/basic/{soil_tag}/`.",
        md,
        count=1,
    )
    # Also retitle the markdown subtitle (van Genuchten–Mualem loam ...) -> soil_tag
    md = md.replace(
        "(van Genuchten\u2013Mualem loam, free-drainage base, single column).",
        f"(van Genuchten\u2013Mualem {soil_tag}, free-drainage base, single column).",
    )
    # Update wet-phase row to reflect compensation (P\u2032 absorbs PE+PT*0.8; PE=PT=0).
    md = re.sub(
        r"\| Wet \(days 5\u201335\) \| P[^|]*\|",
        f"| Wet (days 5\u201335) | P = {forcing['P_WET']:.3f} cm d\u207b\u00b9, PE = {forcing['PE_WET']:.3f}, PT = {forcing['PT_WET']:.3f} cm d\u207b\u00b9 |",
        md,
        count=1,
    )
    # Re-point reference to the renamed intensive-companion notebook.
    md = re.sub(
        r"See `comparison-intensive(?:-[a-z\-]+)?\.ipynb`",
        f"See `comparison-intensive-{soil_tag}.ipynb`",
        md,
    )
    # Strip any prior wet-ET-compensation note (compensation has been retracted).
    md = re.sub(r"\n*\*\*Wet-phase ET compensation:\*\*[\s\S]*?\)\.\n+", "\n\n", md)
    _set_src(cells[0], md)

    # Cell 2: soil + forcing config — replace whole cell.
    # All variant-fixed quantities are loaded from experiment_definitions.json
    # (the single source of truth). Only OAT-tuned MODEL_PARAMS / ET_STRESS
    # remain literal so apply_oat_optima can rewrite them in-place.
    _set_src(cells[2], _make_cell2_basic(soil_tag))

    # Cell 3: HYDRUS-1D setup
    h = _src(cells[3])
    # Output path
    h = h.replace(
        'HYDRUS_WS = HERE / "outputs" / "basic" / "loam" / "phydrus" / "hydrus1d"',
        f'HYDRUS_WS = HERE / "outputs" / "basic" / "{soil_tag}" / "phydrus" / "hydrus1d"',
    )
    h = h.replace(
        'HYDRUS_WS = HERE / "outputs" / "basic" / "phydrus" / "hydrus1d"',
        f'HYDRUS_WS = HERE / "outputs" / "basic" / "{soil_tag}" / "phydrus" / "hydrus1d"',
    )
    # Strip 'loam column experiment' references in comments / descriptions.
    h = h.replace(
        "# Runs HYDRUS-1D for the identical 65-day loam column experiment.",
        f"# Runs HYDRUS-1D for the identical 65-day {soil_tag} column experiment.",
    )
    # Strip the legacy single-material VG_LOAM_CM scalar (above ml = ps.Model(...)).
    h = re.sub(
        r"# ── Convert GWSWEX m-units parameters to cm for HYDRUS-1D[^\n]*\n"
        r"# VG_LOAM format:[^\n]*\n"
        r"VG_LOAM_CM = \[THETA_R, THETA_S, ALPHA / 100\.0, N_VG, K_SAT \* 100\.0, LAM\]\n\n",
        "",
        h,
        count=1,
    )
    # Replace the single-material mat_df registration (below ml.add_waterflow) with a multi-material loop.
    h = re.sub(
        r"mat_df = ml\.get_empty_material_df\(n=1\)\n"
        r"mat_df\.loc\[1\] = VG_LOAM_CM\n"
        r"ml\.add_material\(mat_df\)\n",
        "# Register all soil materials (cm units).\n" + make_hydrus_material_block(),
        h,
        count=1,
    )
    # Replace single-material profile creation with layered profile
    h = re.sub(
        r"# 1-D node profile: 1 cm spacing, hydrostatic IC at WT depth 150 cm\n"
        r"profile = ps\.create_profile\(top=0\.0, bot=-300\.0, dx=1\.0, mat=1\)\n"
        r'profile\["h"\] = -150\.0 - profile\["x"\]\.to_numpy\(float\).*?\n'
        r"ml\.add_profile\(profile\)\n",
        "# 1-D node profile: 1 cm spacing, hydrostatic IC at WT depth (Z_TOP-Z_WT)*100 cm\n"
        "_Z_COL_CM = (Z_TOP - Z_BOT) * 100.0\n"
        "_Z_WT_CM = (Z_TOP - Z_WT) * 100.0\n"
        'profile = ps.create_profile(top=0.0, bot=-_Z_COL_CM, dx=1.0, mat=MATERIALS[0]["id"])\n'
        'profile["h"] = -_Z_WT_CM - profile["x"].to_numpy(float)  # h(z) = -Z_WT_CM - z\n'
        + make_hydrus_profile_per_node_mat_block("Z_TOP * 100.0")
        + "ml.add_profile(profile)\n",
        h,
        count=1,
        flags=re.DOTALL,
    )
    # Description string update
    h = h.replace(
        'description="65-day 3-phase loam column: warmup → wet → dry"',
        f'description="65-day 3-phase {soil_tag} column: warmup → wet → dry"',
    )
    # Feddes poptm must have one entry per material in HYDRUS-1D layered runs.
    h = h.replace(
        "ml.add_root_uptake(model=0, poptm=[-25], p0=-10, p2h=-200, p2l=-800, p3=-8000)",
        "ml.add_root_uptake(model=0, poptm=[-25] * len(MATERIALS), p0=-10, p2h=-200, p2l=-800, p3=-8000)",
    )
    _set_src(cells[3], h)

    # Cell 4: build_model
    b = _src(cells[4])
    b = b.replace(
        'out = HERE / "outputs" / "basic" / "loam" / f"compare-{solver}.nc"',
        f'out = HERE / "outputs" / "basic" / "{soil_tag}" / f"compare-{{solver}}.nc"',
    )
    b = b.replace(
        'out = HERE / "outputs" / "basic" / f"compare-{solver}.nc"',
        f'out = HERE / "outputs" / "basic" / "{soil_tag}" / f"compare-{{solver}}.nc"',
    )
    b = b.replace(
        "m.init_space(ne=1, nl=NL, top=[[Z_TOP]], bot=[list(BNDS[1:])], sID=[[1] * NL], vID=[[1]])",
        "m.init_space(ne=1, nl=NL, top=[[Z_TOP]], bot=[list(BNDS[1:])], sID=[LAYER_MAT_IDS], vID=[[1]])",
    )
    b = b.replace(
        '    m.add_material(id=1, name="loam", K_sat=K_SAT, lam=LAM, vanG=VG)\n',
        "    for _mat in MATERIALS:\n" "        m.add_material(**_mat)\n",
    )
    # Update docstring
    b = b.replace(
        '"""Configure GWSWEXmodel for the 65-day loam column experiment.',
        f'"""Configure GWSWEXmodel for the 65-day {soil_tag} column experiment.',
    )
    _set_src(cells[4], b)

    _apply_post_processing(cells, soil_tag, "basic")

    # Wipe outputs from all code cells (ensure clean variant notebooks)
    for c in cells:
        if c.get("cell_type") == "code":
            c["outputs"] = []
            c["execution_count"] = None

    return nb


# ---------------------------------------------------------------------------
# Intensive notebook transformer
# ---------------------------------------------------------------------------


def transform_intensive(orig: dict, soil_tag: str, layers, forcing: dict) -> dict:
    nb = deepcopy(orig)
    cells = nb["cells"]

    # Cell 0: markdown title
    md = _src(cells[0])
    md = re.sub(
        r"^# .*GWSWEX vs HYDRUS-1D.*$",
        f"# GWSWEX vs HYDRUS-1D: Intensive Comparison — {soil_tag}",
        md,
        count=1,
        flags=re.MULTILINE,
    )
    md = re.sub(
        r"\| Column \| [^|]+\|",
        f"| Column | 1.5 m, 150 × 0.01 m layers, **{soil_tag}** |",
        md,
        count=1,
    )
    md = md.replace(
        "HYDRUS-1D on a 4-phase loam-column experiment.", f"HYDRUS-1D on a 4-phase {soil_tag} column experiment."
    )
    # Wet-phase row: native (uncompensated) hourly forcings.
    _p_cmh = forcing["P_W"] * 100.0  # m h^-1 -> cm h^-1
    _e_cmh = forcing["E_W"] * 100.0
    _t_cmh = forcing["T_W"] * 100.0
    md = re.sub(
        r"\| Wet *\| 10 d *\| [^|]*\| [^|]*\| [^|]*\|",
        f"| Wet      | 10 d     | P = {_p_cmh:.3g} cm h\u207b\u00b9 | E = {_e_cmh:.3g} cm h\u207b\u00b9 | T = {_t_cmh:.3g} cm h\u207b\u00b9 |",
        md,
        count=1,
    )
    # Re-point output path mention to the per-soil subdirectory.
    md = re.sub(
        r"`outputs/intensive(?:/[a-z\-]+)?/`",
        f"`outputs/intensive/{soil_tag}/`",
        md,
    )
    # Strip any prior wet-ET-compensation note (intensive notebooks no longer
    # apply compensation; see docs/sensitivity-analysis-from-comparison.md).
    md = re.sub(r"\n*\*\*Wet-phase ET compensation:\*\*[\s\S]*?\)\.\n?", "", md)
    _set_src(cells[0], md)

    # Cell 2 (index 2): soil + forcing config — replace whole cell.
    # Loaded from experiment_definitions.json; only OAT-tuned MODEL_PARAMS /
    # ET_STRESS remain literal.
    _set_src(cells[2], _make_cell2_intensive(soil_tag))

    # Cell 3 (HYDRUS-1D setup)
    h = _src(cells[3])
    h = h.replace(
        'HYDRUS_WS = HERE / "outputs" / "intensive" / "loam" / "phydrus" / "hydrus1d"',
        f'HYDRUS_WS = HERE / "outputs" / "intensive" / "{soil_tag}" / "phydrus" / "hydrus1d"',
    )
    h = h.replace(
        'HYDRUS_WS = HERE / "outputs" / "intensive" / "phydrus" / "hydrus1d"',
        f'HYDRUS_WS = HERE / "outputs" / "intensive" / "{soil_tag}" / "phydrus" / "hydrus1d"',
    )
    h = h.replace(
        "# Runs HYDRUS-1D for the identical 32-day loam column experiment.",
        f"# Runs HYDRUS-1D for the identical 32-day {soil_tag} column experiment.",
    )
    # Strip legacy single-material VG_LOAM_CM scalars (above ml = ps.Model(...)),
    # keep the geometry conversions Z_COL_CM / Z_WT_CM since later code uses them.
    h = re.sub(
        r"# ── Convert GWSWEX m-units parameters to cm for HYDRUS-1D[^\n]*\n"
        r"# VG_LOAM format:[^\n]*\n"
        r"K_SAT_CMD = K_SAT_MD \* 100\.0[^\n]*\n"
        r"VG_LOAM_CM = \[THETA_R, THETA_S, ALPHA / 100\.0, N_VG, K_SAT_CMD, LAM\]\n\n",
        "",
        h,
        count=1,
    )
    # Replace the single-material mat_df registration with a multi-material loop.
    h = re.sub(
        r"mat_df = ml\.get_empty_material_df\(n=1\)\n"
        r"mat_df\.loc\[1\] = VG_LOAM_CM\n"
        r"ml\.add_material\(mat_df\)\n",
        "# Register all soil materials (cm units).\n" + make_hydrus_material_block(),
        h,
        count=1,
    )
    # Add per-node Mat after the existing profile creation (which already uses Z_COL_CM and Z_WT_CM)
    h = re.sub(
        r"(profile = ps\.create_profile\(top=0\.0, bot=-Z_COL_CM, dx=1\.0, mat=)1\)",
        r'\1MATERIALS[0]["id"])',
        h,
    )
    # Insert per-node Mat assignment after the profile["h"] = ... line
    h = re.sub(
        r"(profile\[\"h\"\] = -Z_WT_CM - profile\[\"x\"\]\.to_numpy\(float\).*?\n)",
        r"\1" + make_hydrus_profile_per_node_mat_block("Z_TOP * 100.0"),
        h,
        count=1,
    )
    # Description
    h = h.replace(
        'description="32-day 4-phase loam column: warmup → wet → dry → cooldown"',
        f'description="32-day 4-phase {soil_tag} column: warmup → wet → dry → cooldown"',
    )
    h = h.replace(
        "ml.add_root_uptake(model=0, poptm=[-25], p0=-10, p2h=-200, p2l=-800, p3=-8000)",
        "ml.add_root_uptake(model=0, poptm=[-25] * len(MATERIALS), p0=-10, p2h=-200, p2l=-800, p3=-8000)",
    )
    _set_src(cells[3], h)

    # Cell 4: build_model
    b = _src(cells[4])
    b = b.replace(
        'out = HERE / "outputs" / "intensive" / "loam" / f"compare-{solver}.nc"',
        f'out = HERE / "outputs" / "intensive" / "{soil_tag}" / f"compare-{{solver}}.nc"',
    )
    b = b.replace(
        'out = HERE / "outputs" / "intensive" / f"compare-{solver}.nc"',
        f'out = HERE / "outputs" / "intensive" / "{soil_tag}" / f"compare-{{solver}}.nc"',
    )
    b = b.replace(
        "m.init_space(ne=NE, nl=NL, top=[[Z_TOP]], bot=bot, sID=[[1] * NL], vID=[[1]])",
        "m.init_space(ne=NE, nl=NL, top=[[Z_TOP]], bot=bot, sID=[LAYER_MAT_IDS], vID=[[1]])",
    )
    b = b.replace(
        '    m.add_material(id=1, name="loam", K_sat=K_SAT_MH, lam=LAM, vanG=VG)\n',
        "    for _mat in MATERIALS:\n"
        '        _mat_h = dict(_mat); _mat_h["K_sat"] = _mat["K_sat"] / 24.0  # m d^-1 -> m h^-1\n'
        "        m.add_material(**_mat_h)\n",
    )
    b = b.replace(
        '"""Configure GWSWEXmodel for the 32-day loam column experiment.',
        f'"""Configure GWSWEXmodel for the 32-day {soil_tag} column experiment.',
    )
    _set_src(cells[4], b)

    _apply_post_processing(cells, soil_tag, "intensive")

    # Wipe outputs from all code cells
    for c in cells:
        if c.get("cell_type") == "code":
            c["outputs"] = []
            c["execution_count"] = None

    return nb


_MP_RE = re.compile(r"MODEL_PARAMS = dict\([^)]*\)")
_ET_RE = re.compile(r"ET_STRESS = dict\([^)]*\)")
_SET_EXPL_RE = re.compile(r'm\.set_solver\(solver="explicit"[^)]*\)')
_SET_IMPL_RE = re.compile(r'm\.set_solver\(solver="implicit"[^)]*\)')


def _fmt_kw(d: dict) -> str:
    parts = []
    for k, v in d.items():
        if isinstance(v, float):
            if v == 0 or 1e-3 <= abs(v) < 1e4:
                parts.append(f"{k}={v:g}")
            else:
                parts.append(f"{k}={v:.0e}")
        else:
            parts.append(f"{k}={v!r}")
    return ", ".join(parts)


def _apply_tuned_params(nb: dict, soil_tag: str, setup: str, tuned: dict) -> None:
    """Inject per-case tuned MODEL_PARAMS + set_solver kwargs in-place."""
    cfg = tuned.get(soil_tag, {}).get(setup, {})
    if not cfg:
        return
    mp_expl = cfg.get("explicit", {}).get("model_params")
    mp_impl = cfg.get("implicit", {}).get("model_params")
    sp_expl = cfg.get("explicit", {}).get("solver_params")
    sp_impl = cfg.get("implicit", {}).get("solver_params")
    # ET stress is per-vegetation; we share a single ET_STRESS dict in the
    # notebook between both solvers. Prefer explicit picks where both differ
    # (explicit is more sensitive to root-zone storage transients via the
    # cascade; implicit smooths via the Picard iteration).
    et_expl = cfg.get("explicit", {}).get("et_stress")
    et_impl = cfg.get("implicit", {}).get("et_stress")
    expl_picks = {p[0] for p in cfg.get("explicit", {}).get("accepted", [])}
    et_merged = dict(et_impl or et_expl or {})
    if et_expl:
        for k in expl_picks & set(et_expl.keys()):
            et_merged[k] = et_expl[k]
    # MODEL_PARAMS is a single global dict shared by both solvers. Merge:
    # if explicit picked a non-default model param value, prefer that (the
    # explicit cascade is the only place ICratio_min is active; implicit is
    # insensitive to it). Otherwise fall back to implicit picks.
    expl_picks = {p[0] for p in cfg.get("explicit", {}).get("accepted", [])}
    mp = dict(mp_impl or mp_expl or {})
    if mp_expl:
        for k in expl_picks & set(mp_expl.keys()):
            mp[k] = mp_expl[k]
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        s = _src(cell)
        s2 = s
        if mp:
            s2 = _MP_RE.sub(f"MODEL_PARAMS = dict({_fmt_kw(mp)})", s2)
        if et_merged:
            s2 = _ET_RE.sub(f"ET_STRESS = dict({_fmt_kw(et_merged)})", s2)
        if sp_expl:
            kw = _fmt_kw(sp_expl)
            s2 = _SET_EXPL_RE.sub(f'm.set_solver(solver="explicit", {kw})', s2)
        if sp_impl:
            kw = _fmt_kw(sp_impl)
            s2 = _SET_IMPL_RE.sub(f'm.set_solver(solver="implicit", {kw})', s2)
        if s2 != s:
            _set_src(cell, s2)


def _compensate_wet_et(fb: dict, fi: dict) -> tuple[dict, dict]:
    """Pass-through: wet-phase ET compensation has been REVERTED everywhere.

    Earlier iterations folded wet-phase ET into precipitation as
    P' = P + 0.8*(PE+PT) with PE_wet = PT_wet = 0 to mask the GWSWEX-vs-
    HYDRUS-1D pre-canopy evaporation accounting mismatch. That artifice
    skewed every metric and damaged the intensive WT trajectory, so it has
    been removed (decision logged in tasks.md P8).
    """
    return dict(fb), dict(fi)


def main() -> None:
    orig_basic = json.loads(SRC_BASIC.read_text())
    orig_int = json.loads(SRC_INT.read_text())

    tuned_path = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "tuned_params.json"
    tuned = json.loads(tuned_path.read_text()) if tuned_path.exists() else {}
    if tuned:
        print(
            f"Loaded tuned params for {sum(len(v) for v in tuned.values())} (soil,setup) pairs from {tuned_path.name}"
        )

    written = []
    for soil_tag, (layers, fb, fi) in VARIANTS.items():
        fb, fi = _compensate_wet_et(fb, fi)
        out_b = EX_DIR / f"comparison-basic-{soil_tag}.ipynb"
        out_i = EX_DIR / f"comparison-intensive-{soil_tag}.ipynb"
        nb_b = transform_basic(orig_basic, soil_tag, layers, fb)
        nb_i = transform_intensive(orig_int, soil_tag, layers, fi)
        _apply_tuned_params(nb_b, soil_tag, "basic", tuned)
        _apply_tuned_params(nb_i, soil_tag, "intensive", tuned)
        out_b.write_text(json.dumps(nb_b, indent=1) + "\n")
        out_i.write_text(json.dumps(nb_i, indent=1) + "\n")
        written.append((out_b.name, out_i.name))

    print(f"Wrote {2 * len(written)} notebooks under {EX_DIR}:")
    for b, i in written:
        print(f"  {b}    {i}")


if __name__ == "__main__":
    main()

"""Update manuscript RMSE/NSE tables and inline summary numbers from validation_summary.md.

Reads:   oat_results/validation_summary.md   (output of validate_all.py)
Writes:  docs/manuscript/manuscript-gwswex.md  (in-place; section tables and summary numbers)

What is updated:
  1. \\label{tab:verif-rmse}  — full-horizon GWH RMSE/NSE (12 rows).
  2. \\label{tab:verif-phase} — phase-resolved GWH RMSE/NSE (24 rows).
  3. Inline summary numbers in the Results text (phase-averaged RMSE/NSE,
     basic-clay-explicit tuned value, etc.).

What is NOT touched:
  - Manuscript structure, prose, references, figures.
  - Any numbers not directly produced by the comparison notebooks.

Usage:
    source .env.d/dev.env && $PYTHON examples/gwswex-vs-hydrus1d/sensitivity-analysis/update_manuscript.py
"""

from __future__ import annotations
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
VAL_SUMMARY = (
    ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "validation_summary.md"
)
TUNED_PARAMS = ROOT / "examples" / "gwswex-vs-hydrus1d" / "sensitivity-analysis" / "oat_results" / "tuned_params.json"
MANUSCRIPT = ROOT / "docs" / "manuscript" / "manuscript-gwswex.md"

# ──────────────────────────────────────────────────────────────────────────────
# Parse validation_summary.md
# ──────────────────────────────────────────────────────────────────────────────

SOIL_DISPLAY = {
    "loam": "loam",
    "sand": "sand",
    "clay": "clay",
    "sand-loam": "sand-over-loam",
    "sand-clay": "sand-over-clay",
    "loam-clay": "loam-over-clay",
}

# Canonical soil order: single-material first (loam, sand, clay), then layered.
# Used consistently across both tables so the rows align.
SOILS_ORDER = ["loam", "sand", "clay", "sand-loam", "sand-clay", "loam-clay"]
SETUPS_ORDER = [("basic", SOILS_ORDER), ("intensive", SOILS_ORDER)]


def _parse_full_table(text: str) -> dict:
    """Return {(setup, soil): {impl/expl: {rmse, nse}}} from the full-horizon table."""
    rows: dict = {}
    in_table = False
    for ln in text.splitlines():
        if "## Full-horizon RMSE and NSE" in ln:
            in_table = True
            continue
        if in_table and ln.startswith("## "):
            break
        if not in_table or not ln.startswith("|"):
            continue
        parts = [p.strip() for p in ln.strip("|").split("|")]
        if len(parts) < 8 or parts[0] in ("setup", "---", ""):
            continue
        setup, soil, rc, dt = parts[0], parts[1], parts[2], parts[3]
        try:
            impl_rmse = float(parts[4]) if parts[4] not in ("FAIL", "-") else None
            impl_nse = float(parts[5]) if parts[5] not in ("FAIL", "-") else None
            expl_rmse = float(parts[6]) if parts[6] not in ("FAIL", "-") else None
            expl_nse = float(parts[7]) if parts[7] not in ("FAIL", "-") else None
        except (ValueError, IndexError):
            continue
        rows[(setup, soil)] = dict(
            impl_rmse=impl_rmse,
            impl_nse=impl_nse,
            expl_rmse=expl_rmse,
            expl_nse=expl_nse,
        )
    return rows


def _parse_phase_table(text: str) -> dict:
    """Return {(setup, soil, solver): {wet_rmse, wet_nse, dry_rmse, dry_nse}}."""
    rows: dict = {}
    in_table = False
    for ln in text.splitlines():
        if "## Phase-resolved RMSE and NSE" in ln:
            in_table = True
            continue
        if in_table and ln.startswith("## "):
            break
        if not in_table or not ln.startswith("|"):
            continue
        parts = [p.strip() for p in ln.strip("|").split("|")]
        if len(parts) < 7 or parts[0] in ("setup", "---", ""):
            continue
        setup, soil, slv = parts[0], parts[1], parts[2]
        # solver codes in the MD are "impl" / "expl" — map to manuscript terms
        solver = "implicit" if slv == "impl" else "explicit"
        try:
            wet_rmse = float(parts[3]) if parts[3] not in ("FAIL", "-") else None
            wet_nse = float(parts[4]) if parts[4] not in ("FAIL", "-") else None
            dry_rmse = float(parts[5]) if parts[5] not in ("FAIL", "-") else None
            dry_nse = float(parts[6]) if parts[6] not in ("FAIL", "-") else None
        except (ValueError, IndexError):
            continue
        rows[(setup, soil, solver)] = dict(
            wet_rmse=wet_rmse,
            wet_nse=wet_nse,
            dry_rmse=dry_rmse,
            dry_nse=dry_nse,
        )
    return rows


# ──────────────────────────────────────────────────────────────────────────────
# LaTeX table builders
# ──────────────────────────────────────────────────────────────────────────────

BASIC_SOILS = SOILS_ORDER
INT_SOILS = SOILS_ORDER

_DASH = "--"  # LaTeX em-dash placeholder for missing NSE


def _nse(v) -> str:
    if v is None:
        return _DASH
    return f"{v:.3f}"


def _rmse(v) -> str:
    if v is None:
        return _DASH
    return f"{v:.2f}"


def build_verif_rmse_table(full: dict) -> str:
    """Produce the tabular body rows for tab:verif-rmse."""
    lines = []
    for setup, soils in SETUPS_ORDER:
        for soil in soils:
            d = full.get((setup, soil), {})
            disp = SOIL_DISPLAY.get(soil, soil)
            lines.append(
                f"{setup} & {disp} & "
                f"{_rmse(d.get('impl_rmse'))} & {_nse(d.get('impl_nse'))} & "
                f"{_rmse(d.get('expl_rmse'))} & {_nse(d.get('expl_nse'))} \\\\"
            )
    return "\n".join(lines)


def build_verif_phase_table(phase: dict) -> str:
    """Produce the tabular body rows for tab:verif-phase."""
    lines = []
    for setup, soils in SETUPS_ORDER:
        for soil in soils:
            disp = SOIL_DISPLAY.get(soil, soil)
            for solver in ("implicit", "explicit"):
                slv_abbr = "impl." if solver == "implicit" else "expl."
                d = phase.get((setup, soil, solver), {})
                lines.append(
                    f"{setup}, {disp} & {slv_abbr} & "
                    f"{_rmse(d.get('wet_rmse'))} & {_nse(d.get('wet_nse'))} & "
                    f"{_rmse(d.get('dry_rmse'))} & {_nse(d.get('dry_nse'))} & "
                    f"{_rmse(None)} & {_nse(None)} \\\\"  # total columns filled from full table
                )
    return "\n".join(lines)


def build_verif_phase_table_v2(phase: dict, full: dict) -> str:
    """Same as build_verif_phase_table but also fills total RMSE/NSE from full dict."""
    lines = []
    for setup, soils in SETUPS_ORDER:
        for soil in soils:
            disp = SOIL_DISPLAY.get(soil, soil)
            fd = full.get((setup, soil), {})
            for solver in ("implicit", "explicit"):
                slv_abbr = "impl." if solver == "implicit" else "expl."
                d = phase.get((setup, soil, solver), {})
                if solver == "implicit":
                    tot_rmse = fd.get("impl_rmse")
                    tot_nse = fd.get("impl_nse")
                else:
                    tot_rmse = fd.get("expl_rmse")
                    tot_nse = fd.get("expl_nse")
                lines.append(
                    f"{setup}, {disp} & {slv_abbr} & "
                    f"{_rmse(d.get('wet_rmse'))} & {_nse(d.get('wet_nse'))} & "
                    f"{_rmse(d.get('dry_rmse'))} & {_nse(d.get('dry_nse'))} & "
                    f"{_rmse(tot_rmse)} & {_nse(tot_nse)} \\\\"
                )
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────────────
# Inline-number updaters
# ──────────────────────────────────────────────────────────────────────────────


def _mean(values) -> float | None:
    vals = [v for v in values if v is not None]
    return sum(vals) / len(vals) if vals else None


def compute_phase_averages(phase: dict, full: dict) -> dict[str, float | None]:
    """Compute the summary statistics cited in the Results text."""
    basic_impl_wet = [phase.get(("basic", s, "implicit"), {}).get("wet_rmse") for s in BASIC_SOILS]
    basic_expl_wet = [phase.get(("basic", s, "explicit"), {}).get("wet_rmse") for s in BASIC_SOILS]
    basic_impl_dry = [phase.get(("basic", s, "implicit"), {}).get("dry_rmse") for s in BASIC_SOILS]
    basic_expl_dry = [phase.get(("basic", s, "explicit"), {}).get("dry_rmse") for s in BASIC_SOILS]
    int_impl_wet = [phase.get(("intensive", s, "implicit"), {}).get("wet_rmse") for s in INT_SOILS]
    int_expl_wet = [phase.get(("intensive", s, "explicit"), {}).get("wet_rmse") for s in INT_SOILS]
    int_impl_dry = [phase.get(("intensive", s, "implicit"), {}).get("dry_rmse") for s in INT_SOILS]
    int_expl_dry = [phase.get(("intensive", s, "explicit"), {}).get("dry_rmse") for s in INT_SOILS]

    # Look up the OAT-run baseline RMSE for basic-clay-explicit from tuned_params.json
    # so the inline sentence "from X cm at the notebook baseline" stays accurate
    # after the forcing change.
    clay_expl_baseline: float | None = None
    if TUNED_PARAMS.exists():
        try:
            tp = json.loads(TUNED_PARAMS.read_text())
            clay_expl_baseline = tp.get("clay", {}).get("basic", {}).get("explicit", {}).get("baseline_rmse")
        except Exception:
            pass

    return dict(
        basic_impl_wet=_mean(basic_impl_wet),
        basic_expl_wet=_mean(basic_expl_wet),
        basic_impl_dry=_mean(basic_impl_dry),
        basic_expl_dry=_mean(basic_expl_dry),
        int_impl_wet=_mean(int_impl_wet),
        int_expl_wet=_mean(int_expl_wet),
        int_impl_dry=_mean(int_impl_dry),
        int_expl_dry=_mean(int_expl_dry),
        # Individual cells for inline mentions
        clay_basic_expl_full=full.get(("basic", "clay"), {}).get("expl_rmse"),
        clay_basic_expl_baseline=clay_expl_baseline,
    )


_BASIC_AVG_RE = re.compile(
    r"Averaged across the six basic-setup soils, the wet-phase GWH-RMSE is "
    r"[\d.]+ cm \(implicit\) and [\d.]+ cm \(explicit\) and the dry-phase "
    r"GWH-RMSE is [\d.]+ cm \(implicit\) and [\d.]+ cm \(explicit\)\."
)

_INT_AVG_RE = re.compile(
    r"averaged across soils the wet-phase RMSE is [\d.]+ cm \(implicit\) and "
    r"[\d.]+ cm \(explicit\), against a dry-phase RMSE of [\d.]+ cm \(implicit\) "
    r"and [\d.]+ cm \(explicit\)"
)

_CLAY_EXPL_RE = re.compile(
    r"lifts the GWH RMSE from approximately [\d.]+ cm at the notebook baseline "
    r"to ([\d.]+) cm in the OAT-tuned configuration"
)


def update_inline_numbers(text: str, avgs: dict, full: dict) -> str:
    """Replace computed inline summary numbers in the manuscript text."""
    biw = avgs.get("basic_impl_wet")
    bew = avgs.get("basic_expl_wet")
    bid = avgs.get("basic_impl_dry")
    bed = avgs.get("basic_expl_dry")
    iiw = avgs.get("int_impl_wet")
    iew = avgs.get("int_expl_wet")
    iid = avgs.get("int_impl_dry")
    ied = avgs.get("int_expl_dry")
    clay_expl = avgs.get("clay_basic_expl_full")

    if all(v is not None for v in [biw, bew, bid, bed]):
        replacement = (
            f"Averaged across the six basic-setup soils, the wet-phase GWH-RMSE is "
            f"{biw:.2f} cm (implicit) and {bew:.2f} cm (explicit) and the dry-phase "
            f"GWH-RMSE is {bid:.2f} cm (implicit) and {bed:.2f} cm (explicit)."
        )
        text = _BASIC_AVG_RE.sub(replacement, text)
        print(f"  Updated basic-setup average sentence.")

    if all(v is not None for v in [iiw, iew, iid, ied]):
        replacement = (
            f"averaged across soils the wet-phase RMSE is {iiw:.2f} cm (implicit) and "
            f"{iew:.2f} cm (explicit), against a dry-phase RMSE of {iid:.2f} cm (implicit) "
            f"and {ied:.2f} cm (explicit)"
        )
        text = _INT_AVG_RE.sub(replacement, text)
        print(f"  Updated intensive-setup average sentence.")

    if clay_expl is not None:
        clay_base = avgs.get("clay_basic_expl_baseline")
        base_str = f"{clay_base:.1f}" if clay_base is not None else "11.9"
        text = _CLAY_EXPL_RE.sub(
            lambda m: f"lifts the GWH RMSE from approximately {base_str} cm at the notebook baseline "
            f"to {clay_expl:.2f} cm in the OAT-tuned configuration",
            text,
        )
        print(f"  Updated basic-clay-explicit: baseline={base_str} cm → tuned={clay_expl:.2f} cm.")

    return text


_INT_CLAY_CEIL_RE = re.compile(r"The intensive-clay case \(both solvers, [\d]+-+[\d]+ cm GWH RMSE\)")
_INT_LOAM_CEIL_RE = re.compile(
    r"The intensive-loam and intensive-sand-over-loam cases \(both solvers, [\d]+-+[\d]+ cm GWH RMSE\)"
)


def _range_str(values: list) -> str:
    """Format a list of RMSE floats as an integer 'min--max' LaTeX range string."""
    import math

    valid = [v for v in values if v is not None]
    if not valid:
        return "?--?"
    lo = math.floor(min(valid))
    hi = math.ceil(max(valid))
    return f"{lo}--{hi}"


def update_structural_ceilings(text: str, full: dict) -> str:
    """Update the 'X--Y cm GWH RMSE' structural ceiling ranges in the manuscript."""
    clay_vals = [
        full.get(("intensive", "clay"), {}).get("impl_rmse"),
        full.get(("intensive", "clay"), {}).get("expl_rmse"),
    ]
    loam_sl_vals = [
        full.get(("intensive", "loam"), {}).get("impl_rmse"),
        full.get(("intensive", "loam"), {}).get("expl_rmse"),
        full.get(("intensive", "sand-loam"), {}).get("impl_rmse"),
        full.get(("intensive", "sand-loam"), {}).get("expl_rmse"),
    ]

    clay_range = _range_str(clay_vals)
    loam_sl_range = _range_str(loam_sl_vals)

    text, n1 = _INT_CLAY_CEIL_RE.subn(f"The intensive-clay case (both solvers, {clay_range} cm GWH RMSE)", text)
    text, n2 = _INT_LOAM_CEIL_RE.subn(
        f"The intensive-loam and intensive-sand-over-loam cases (both solvers, {loam_sl_range} cm GWH RMSE)", text
    )
    if n1:
        print(f"  Updated intensive-clay structural ceiling: {clay_range} cm.")
    else:
        print("  WARNING: intensive-clay ceiling pattern not found.")
    if n2:
        print(f"  Updated intensive-loam/sand-loam structural ceiling: {loam_sl_range} cm.")
    else:
        print("  WARNING: intensive-loam/sand-loam ceiling pattern not found.")

    return text


# ──────────────────────────────────────────────────────────────────────────────
# Table replacement in manuscript
# ──────────────────────────────────────────────────────────────────────────────

_TABLE_ROW_RE = re.compile(
    r"(basic|intensive) & (loam|sand|clay|sand-over-loam|sand-over-clay|loam-over-clay)"
    r" & [\d.]+|[\-\-]+ & [\d.]+|[\-\-]+ & [\d.]+|[\-\-]+ & [\d.]+|[\-\-]+"
)


def _replace_tabular_body(text: str, label: str, new_body: str) -> str:
    """Replace the row content inside a table that ends with \\label{<label>}."""
    # Find the tabular block containing the label.
    # Strategy: locate \hline ... \hline block preceding \end{tabular} for this table.
    # The table label appears in the \caption{} line.
    label_pat = re.compile(r"\\label\{" + re.escape(label) + r"\}")
    m_label = label_pat.search(text)
    if not m_label:
        print(f"  WARNING: label {label} not found in manuscript. Skipping.")
        return text

    # In the manuscript, \label{} sits inside \caption{} which appears *before*
    # the \begin{tabular} block.  Both forward-searches must start after the label.
    tab_start = text.find(r"\begin{tabular}", m_label.start())
    tab_end = text.find(r"\end{tabular}", m_label.start())
    if tab_start == -1 or tab_end == -1:
        print(f"  WARNING: could not bound tabular block for {label}. Skipping.")
        return text

    block = text[tab_start : tab_end + len(r"\end{tabular}")]

    # Split block on the second \hline (end of header row).
    hline_positions = [m.start() for m in re.finditer(r"\\hline", block)]
    if len(hline_positions) < 2:
        print(f"  WARNING: tabular block for {label} has fewer than 2 \\hlines. Skipping.")
        return text

    # Row body starts after the second \hline (position hline_positions[1] + len('\hline')).
    header_end = hline_positions[1] + len(r"\hline")
    # Row body ends at the last \hline (before \end{tabular}).
    last_hline = hline_positions[-1]

    old_body = block[header_end:last_hline]
    new_block = block[:header_end] + "\n" + new_body + "\n" + block[last_hline:]
    print(f"  Replaced {label} tabular body ({len(old_body)} → {len(new_block) - len(block) + len(old_body)} chars).")
    return text[:tab_start] + new_block + text[tab_end + len(r"\end{tabular}") :]


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────


def main() -> None:
    if not VAL_SUMMARY.exists():
        raise SystemExit(f"validation_summary.md not found: {VAL_SUMMARY}\nRun validate_all.py first.")

    summary_text = VAL_SUMMARY.read_text()
    full = _parse_full_table(summary_text)
    phase = _parse_phase_table(summary_text)
    print(f"Parsed {len(full)} full-horizon and {len(phase)} phase-resolved entries.")

    manuscript = MANUSCRIPT.read_text()

    # 1. Replace tab:verif-rmse rows.
    new_verif_rows = build_verif_rmse_table(full)
    manuscript = _replace_tabular_body(manuscript, "tab:verif-rmse", new_verif_rows)

    # 2. Replace tab:verif-phase rows.
    new_phase_rows = build_verif_phase_table_v2(phase, full)
    manuscript = _replace_tabular_body(manuscript, "tab:verif-phase", new_phase_rows)

    # 3. Update inline summary numbers.
    avgs = compute_phase_averages(phase, full)
    manuscript = update_inline_numbers(manuscript, avgs, full)

    # 4. Update structural ceiling ranges.
    manuscript = update_structural_ceilings(manuscript, full)

    MANUSCRIPT.write_text(manuscript)
    print(f"\nUpdated {MANUSCRIPT.name}.")


if __name__ == "__main__":
    main()

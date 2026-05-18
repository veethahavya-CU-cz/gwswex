"""Update docs/sensitivity-analysis-from-comparison.md from post-OAT pipeline outputs.

Reads:   oat_results/validation_summary.md   (output of validate_all.py)
         oat_results/tuned_params.json        (output of apply_oat_optima.py)
Writes:  docs/sensitivity-analysis-from-comparison.md  (in-place)

What is updated:
  - §3.1 Table 3.1: Full-horizon RMSE/NSE (12 rows).
  - §3.2 Implicit solver table: Final OAT-tuned implicit parameters.
  - §3.2 Explicit solver table: Final OAT-tuned explicit parameters.
  - §4.3 inline references: "11.9 cm to 2.36 cm" and other per-cell RMSE numbers.

What is NOT touched:
  - Prose sections §1, §2, §4.1, §4.2 (structural analysis — not results-dependent).
  - §4.3 qualitative mechanism text (tuning-direction sentences are results-independent).
  - References, headings, or any other structural element.

Usage:
    source .env.d/dev.env
    $PYTHON examples/gwswex-vs-hydrus1d/sensitivity-analysis/update_companion_doc.py
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
COMPANION = ROOT / "docs" / "sensitivity-analysis-from-comparison.md"

SOILS_ORDER = ["loam", "sand", "clay", "sand-loam", "sand-clay", "loam-clay"]
SETUPS_ORDER = ["basic", "intensive"]


# ──────────────────────────────────────────────────────────────────────────────
# Parse validation_summary.md (full-horizon table only)
# ──────────────────────────────────────────────────────────────────────────────


def _parse_full_table(text: str) -> dict:
    """Return {(setup, soil): {impl_rmse, impl_nse, expl_rmse, expl_nse}}."""
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
        setup, soil = parts[0], parts[1]
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


# ──────────────────────────────────────────────────────────────────────────────
# §3.1 Table 3.1 replacement (markdown pipe table)
# ──────────────────────────────────────────────────────────────────────────────

_HEADER_31 = "| Setup | Soil | Implicit RMSE [cm] | Implicit NSE | Explicit RMSE [cm] | Explicit NSE |"
_SEP_31 = "|---|---|---|---|---|---|"


def _nse_str(v) -> str:
    if v is None:
        return "--"
    if v < 0:
        return f"${v:.3f}$"
    return f"{v:.3f}"


def _rmse_str(v) -> str:
    return f"{v:.2f}" if v is not None else "--"


def build_table31(full: dict) -> str:
    lines = [_HEADER_31, _SEP_31]
    for setup in SETUPS_ORDER:
        for soil in SOILS_ORDER:
            d = full.get((setup, soil), {})
            lines.append(
                f"| {setup} | {soil} | "
                f"{_rmse_str(d.get('impl_rmse'))} | {_nse_str(d.get('impl_nse'))} | "
                f"{_rmse_str(d.get('expl_rmse'))} | {_nse_str(d.get('expl_nse'))} |"
            )
    return "\n".join(lines)


def replace_table31(text: str, new_table: str) -> str:
    """Replace Table 3.1 in the companion doc."""
    # The table starts with the static preamble text ending in "reproduced in Table 3.1."
    # then a blank line, then the header row.
    anchor = "**Table 3.1.**"
    m = re.search(re.escape(anchor), text)
    if not m:
        print("  WARNING: Table 3.1 anchor not found. Skipping.")
        return text

    # Find the header line after the anchor
    header_start = text.find(_HEADER_31, m.start())
    if header_start == -1:
        print("  WARNING: Table 3.1 header not found. Skipping.")
        return text

    # Find the end of the table: first blank line after the header
    after_header = text.find("\n", header_start)
    after_sep = text.find("\n", after_header + 1)  # past separator

    # Scan for end of table: first line that doesn't start with '|'
    pos = after_sep + 1
    while pos < len(text):
        eol = text.find("\n", pos)
        if eol == -1:
            eol = len(text)
        line = text[pos:eol].strip()
        if line and not line.startswith("|"):
            break
        pos = eol + 1

    # pos now points to the start of the line after the last table row
    old_table = text[header_start:pos].rstrip()
    print(f"  Replaced Table 3.1 ({len(old_table)} chars).")
    return text[:header_start] + new_table + "\n" + text[pos:]


# ──────────────────────────────────────────────────────────────────────────────
# §3.2 Implicit solver table replacement
# ──────────────────────────────────────────="──
# ──────────────────────────────────────────────────────────────────────────────

_IMPL_HEADER = (
    "| Soil | Setup | RMSE [cm] | psi_f | F_min | ICratio_min | "
    "$s^*$ | $s_w$ | $s_h$ | $s_e$ | picard_tol | picard_max_iter | n_trapz | beta_hyst |"
)
_EXPL_HEADER = (
    "| Soil | Setup | RMSE [cm] | psi_f | F_min | ICratio_min | "
    "$s^*$ | $s_w$ | $s_h$ | $s_e$ | courant_number | n_trapz | beta_hyst |"
)
_IMPL_SEP = "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
_EXPL_SEP = "|---|---|---|---|---|---|---|---|---|---|---|---|---|"

# Sentinel for parameters that don't apply to a solver (e.g. n_trapz, beta_hyst
# are not in all solver param dicts for the intensive implicit solver).
_NA = "--"


def _pval(v) -> str:
    """Format a parameter value for the table."""
    if v is None:
        return _NA
    if isinstance(v, float):
        # Use scientific notation for very small numbers; plain float otherwise.
        if abs(v) < 1e-3 and v != 0.0:
            return f"{v:.0e}"
        return str(v)
    return str(v)


def build_impl_table32(tuned: dict, full: dict) -> str:
    lines = [_IMPL_HEADER, _IMPL_SEP]
    for setup in SETUPS_ORDER:
        for soil in SOILS_ORDER:
            d = tuned.get(soil, {}).get(setup, {}).get("implicit", {})
            fd = full.get((setup, soil), {})
            rmse = _rmse_str(d.get("final_rmse") or fd.get("impl_rmse"))
            sp = d.get("solver_params", {})
            mp = d.get("model_params", {})
            et = d.get("et_stress", {})
            lines.append(
                f"| {soil} | {setup} | {rmse} | "
                f"{_pval(mp.get('psi_f'))} | {_pval(mp.get('F_min'))} | "
                f"{_pval(mp.get('ICratio_min'))} | "
                f"{_pval(et.get('s_star'))} | {_pval(et.get('s_w'))} | "
                f"{_pval(et.get('s_h'))} | {_pval(et.get('s_e'))} | "
                f"{_pval(sp.get('picard_tol'))} | {_pval(sp.get('picard_max_iter'))} | "
                f"{_pval(sp.get('n_trapz', _NA))} | {_pval(sp.get('beta_hyst', _NA))} |"
            )
    return "\n".join(lines)


def build_expl_table32(tuned: dict, full: dict) -> str:
    lines = [_EXPL_HEADER, _EXPL_SEP]
    for setup in SETUPS_ORDER:
        for soil in SOILS_ORDER:
            d = tuned.get(soil, {}).get(setup, {}).get("explicit", {})
            fd = full.get((setup, soil), {})
            rmse = _rmse_str(d.get("final_rmse") or fd.get("expl_rmse"))
            sp = d.get("solver_params", {})
            mp = d.get("model_params", {})
            et = d.get("et_stress", {})
            lines.append(
                f"| {soil} | {setup} | {rmse} | "
                f"{_pval(mp.get('psi_f'))} | {_pval(mp.get('F_min'))} | "
                f"{_pval(mp.get('ICratio_min'))} | "
                f"{_pval(et.get('s_star'))} | {_pval(et.get('s_w'))} | "
                f"{_pval(et.get('s_h'))} | {_pval(et.get('s_e'))} | "
                f"{_pval(sp.get('courant_number'))} | "
                f"{_pval(sp.get('n_trapz', _NA))} | {_pval(sp.get('beta_hyst', _NA))} |"
            )
    return "\n".join(lines)


def _replace_solver_table(text: str, header: str, sep: str, new_body: str, label: str) -> str:
    """Replace the data rows of the implicit or explicit solver table in §3.2."""
    m = re.search(re.escape(header), text)
    if not m:
        print(f"  WARNING: {label} table header not found. Skipping.")
        return text

    # Find the separator row after the header
    sep_start = text.find(sep, m.end())
    if sep_start == -1:
        print(f"  WARNING: {label} table separator not found. Skipping.")
        return text

    # Find the end of the separator row
    after_sep = text.find("\n", sep_start) + 1

    # Scan rows until a non-pipe line
    pos = after_sep
    while pos < len(text):
        eol = text.find("\n", pos)
        if eol == -1:
            eol = len(text)
        line = text[pos:eol].strip()
        if line and not line.startswith("|"):
            break
        pos = eol + 1

    old_rows = text[after_sep:pos].rstrip()
    print(f"  Replaced {label} table ({len(old_rows)} → {len(new_body)} chars).")
    return text[:after_sep] + new_body + "\n" + text[pos:]


# ──────────────────────────────────────────────────────────────────────────────
# §4.3 inline RMSE references
# ──────────────────────────────────────────────────────────────────────────────

_CLAY_EXPL_RE = re.compile(
    r"reducing `courant_number` from \S+ to \S+, paired with `n_trapz = \d+"
    r"`, "
    r"reduces the GWH RMSE from approximately [\d.]+ cm to ([\d.]+) cm"
)


def update_inline_43(text: str, tuned: dict, full: dict) -> str:
    """Update specific RMSE numbers in §4.3 discussion."""
    clay_expl = full.get(("basic", "clay"), {}).get("expl_rmse")
    clay_base = tuned.get("clay", {}).get("basic", {}).get("explicit", {}).get("baseline_rmse")

    if clay_expl is not None:
        base_str = f"{clay_base:.1f}" if clay_base is not None else "11.9"

        def _clay_sub(m):
            # Preserve the courant_number and n_trapz values from the matched text.
            return (
                m.group(0).split("reduces the GWH RMSE")[0]
                + f"reduces the GWH RMSE from approximately {base_str} cm to {clay_expl:.2f} cm"
            )

        text, n = _CLAY_EXPL_RE.subn(_clay_sub, text)
        if n:
            print(f"  Updated §4.3 basic-clay-explicit: baseline={base_str} → tuned={clay_expl:.2f} cm.")
        else:
            print("  WARNING: §4.3 basic-clay-explicit pattern not found.")
    return text


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────


def main() -> None:
    if not VAL_SUMMARY.exists():
        raise SystemExit(f"validation_summary.md not found: {VAL_SUMMARY}\nRun validate_all.py first.")
    if not TUNED_PARAMS.exists():
        raise SystemExit(f"tuned_params.json not found: {TUNED_PARAMS}\nRun apply_oat_optima.py first.")

    summary_text = VAL_SUMMARY.read_text()
    full = _parse_full_table(summary_text)
    tuned = json.loads(TUNED_PARAMS.read_text())
    print(f"Parsed {len(full)} full-horizon validation entries.")
    print(f"Loaded tuned params for {len(tuned)} soils.")

    companion = COMPANION.read_text()

    # 1. §3.1 Table 3.1
    new_t31 = build_table31(full)
    companion = replace_table31(companion, new_t31)

    # 2. §3.2 Implicit solver table
    new_impl = build_impl_table32(tuned, full)
    companion = _replace_solver_table(
        companion,
        _IMPL_HEADER,
        _IMPL_SEP,
        new_impl.split("\n", 2)[2],  # skip header+sep rows (already in doc)
        "implicit solver",
    )

    # 3. §3.2 Explicit solver table
    new_expl = build_expl_table32(tuned, full)
    companion = _replace_solver_table(
        companion, _EXPL_HEADER, _EXPL_SEP, new_expl.split("\n", 2)[2], "explicit solver"  # skip header+sep rows
    )

    # 4. §4.3 inline RMSE
    companion = update_inline_43(companion, tuned, full)

    COMPANION.write_text(companion)
    print(f"\nUpdated {COMPANION.name}.")


if __name__ == "__main__":
    main()

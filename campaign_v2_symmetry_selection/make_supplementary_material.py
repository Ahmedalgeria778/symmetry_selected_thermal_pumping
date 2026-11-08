#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble supplementary_material/ from the data actually used by main_v2.tex.

The manuscript points to "the supplementary material" for the dense Pe_c survey
(173 wavenumbers x 5 Peclet numbers) and for the campaign data behind the
reversal map. This script gathers those files, records a manifest with byte
counts, SHA-256 digests and row counts, and writes the README that describes
each item and the command that regenerates it. Nothing is invented here: every
file is either copied unchanged or written by this script from the copied data.

    py -X utf8 campaign_v2_symmetry_selection/make_supplementary_material.py
"""
from __future__ import annotations

import hashlib
import os
import shutil
import sys
from pathlib import Path

WORK = Path(r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification")
SUP = WORK / "supplementary_material"
OUT = WORK / "campaign_v2_symmetry_selection"
SUITE = WORK / "dns_physical_validation_suite"
SURVEY = WORK / "validation_outputs_ULTIMATE_DENSE_CPU"

GROUPS = [
    ("A_pe_c_dense_survey", "S1",
     "Dense finite-domain Pe_c survey: 173 wavenumbers on 0.02<=k<=10 for five "
     "Pe values (865 second-order evaluations) at Re=10, Gr=100, Pr=0.7. This is "
     "the data behind the endpoint classification quoted in the main text; the "
     "main text retains only the two bracketing calculations Pe=10.85 and "
     "Pe=10.9.",
     [SURVEY / f for f in (
         "ultimate_dense_cpu_scan.csv", "ultimate_Pe_summary.csv",
         "ultimate_roots.csv", "ultimate_sign_check.csv",
         "ultimate_local_boundary_check.csv", "ULTIMATE_CPU_VALIDATION.json",
         "table_ultimate_global_summary.tex", "table_ultimate_roots.tex",
         "table_ultimate_sign_check.tex")]),
    ("B_perturbative_campaigns", "S2",
     "Reversal cartography and phase-selection data: the D.1 Reynolds scan "
     "(primary window Re<=40 and the low-k extension to Re=120), the Brent-refined "
     "crossings, the power-law fits, the D.2 Grashof-independence test and the "
     "D.5/C.2 phase-law outputs.",
     [OUT / f for f in (
         "D1_Re_grid.csv", "D1_Re_grid_lowk.csv", "D1_Re_c.csv",
         "D1_Re_c_all.csv", "D1_Re_c_fit.json",
         "D2_Gr_invariance.csv", "D2_Re_c_vs_Gr.csv", "D2_Re_c_spread.json",
         "D5_phase_law.csv", "D5_phase_sweep_direct.csv",
         "D5_reference_coefficients.json", "D5_ratio_map.csv",
         "P2_Gr_linearity_two_coeffs.csv", "P2_Re_scan_two_coeffs.csv",
         "P2_two_coeff_law_validation.csv", "P2_Re_c_two_branches.json",
         "C2_phase_lag.csv", "C2_phase_unwrapped.csv")]),
    ("C_dns_validation", "S3",
     "Nonlinear DNS results reported in the main text (amplitude, "
     "spatial and time-step convergence, Grashof scaling, fixed-Pe Reynolds "
     "sweep, reflection and symmetric-forcing tests, master summary) and the "
     "completed D.3/D.4/D.5 campaign: pre-run predictions, raw per-case "
     "measurements with residuals and wall times, the derived deviation table "
     "and the analysis summary behind the DNS confirmation section, plus the "
     "Grashof-pair matrix and the reduced data of its figure, which test "
     "Corollary 1 nonlinearly.",
     [SUITE / f for f in (
         "dns_amplitude_convergence.csv", "dns_spatial_convergence.csv",
         "dns_timestep_convergence.csv", "dns_gr_scaling.csv",
         "dns_re_sweep.csv", "dns_symmetry_test.csv",
         "dns_master_validation.csv")]
     + [OUT / f for f in (
         "DNS_D345_predictions.csv", "DNS_D345_results.csv",
         "DNS_D345_measured.csv", "DNS_D345_analysis.json",
         "DNS_D2_grpair_results.csv", "DNS_D2_grpair_measured.csv",
         "DNS_D2_grpair_analysis.json", "FIG_v2_11_dns_grpair.json")]),
    ("D_scripts", "S4",
     "Code that produced every file above. The solver is the CPU-authoritative "
     "SciPy solve_bvp chain; run_campaigns.py generates group S2, "
     "run_D1_lowk_extension.py the low-k extension and the locus figures, "
     "ultimate_dense_cpu_validation_publication_CORRECTED.py group S1, "
     "dns_physical_validation_suite.py group S3, and run_dns_D3_D4_D5.py the "
     "DNS matrix (dry run by default, --go to integrate) with "
     "analyze_dns_D345.py turning its checkpoint into the measured table, the "
     "deviation analysis and the multi-wavenumber figure. "
     "run_dns_D2_grpair.py runs the four-case Grashof-pair matrix and "
     "analyze_dns_D2_grpair.py reduces it to the measured Grashof ratios. "
     "verify_manuscript_numbers.py re-derives every number quoted in the "
     "manuscript, and check_build.py plus check_cover_letter.py gate the two "
     "compiled documents.",
     [WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py",
      WORK / "dns_physical_validation_suite.py",
      WORK / "ultimate_dense_cpu_validation_publication_CORRECTED.py",
      OUT / "run_campaigns.py", OUT / "run_D1_lowk_extension.py",
      OUT / "run_phase_direct_sweep.py", OUT / "run_phase_clean.py",
      OUT / "run_quadrature.py", OUT / "run_dns_D3_D4_D5.py",
      OUT / "run_dns_D2_grpair.py", OUT / "analyze_dns_D2_grpair.py",
      OUT / "analyze_dns_D345.py", OUT / "smoke_dns_import.py",
      OUT / "launch_dns_D345.ps1", OUT / "watch_dns_D345.ps1",
      OUT / "make_fig_endpoint.py", OUT / "verify_manuscript_numbers.py",
      OUT / "check_notation.py", OUT / "check_duplicate_text.py",
      OUT / "check_build.py", OUT / "check_cover_letter.py",
      OUT / "make_fig_dns_grpair.py"]),
]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def rows(path):
    if path.suffix.lower() not in (".csv", ".tex", ".py", ".json", ".md"):
        return ""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh)


def main():
    if SUP.exists():
        shutil.rmtree(SUP)
    manifest = []
    missing = []
    for group, tag, _desc, files in GROUPS:
        dest_dir = SUP / group
        dest_dir.mkdir(parents=True, exist_ok=True)
        for src in files:
            if not src.exists():
                missing.append(str(src))
                continue
            dest = dest_dir / src.name
            shutil.copy2(src, dest)
            manifest.append(dict(group=group, tag=tag, file=src.name,
                                 bytes=dest.stat().st_size,
                                 lines=rows(dest), sha256=digest(dest),
                                 origin=str(src.relative_to(WORK))))
            print("copied %-6s %s" % (tag, src.name))
    import pandas as pd
    mdf = pd.DataFrame(manifest)
    mdf.to_csv(SUP / "MANIFEST.csv", index=False)
    copied = set(mdf.file)

    readme = ["# Supplementary material - Symmetry-selected nonlinear thermal "
              "pumping (main_v2.tex)", "",
              "Assembled by `campaign_v2_symmetry_selection/"
              "make_supplementary_material.py`; every file is copied unchanged "
              "from the project tree and listed with its byte count, line count "
              "and SHA-256 digest in `MANIFEST.csv`. No file here contains data "
              "that is not used or not quoted by the manuscript.", ""]
    for group, tag, desc, files in GROUPS:
        readme += ["## %s - %s" % (tag, group.replace("_", " ")), "", desc, ""]
        for src in files:
            if src.name in copied:
                readme.append("* `%s`" % src.name)
        readme.append("")
    readme += ["## Reproducing the manuscript numbers from this package", "",
               "The scripts are provided for provenance and are the exact files "
               "used in the project. They address the project tree through "
               "absolute paths, so regeneration means running them from the "
               "repository root, not from inside this package:", "",
               "```",
               "py -X utf8 campaign_v2_symmetry_selection/verify_manuscript_numbers.py",
               "py -X utf8 campaign_v2_symmetry_selection/run_D1_lowk_extension.py",
               "py -X utf8 campaign_v2_symmetry_selection/run_campaigns.py",
               "py -X utf8 campaign_v2_symmetry_selection/run_dns_D3_D4_D5.py   # dry run only",
               "py -X utf8 campaign_v2_symmetry_selection/check_notation.py",
               "```", "",
               "The pending DNS confirmations (D.3/D.4/D.5) are described in the "
               "Outlook subsection of the manuscript; their predictions are in "
               "`C_dns_validation/DNS_D345_predictions.csv` and are regenerated "
               "by the dry run of `D_scripts/run_dns_D3_D4_D5.py`."]
    (SUP / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print("\n%d files, %s README + MANIFEST" % (len(manifest), SUP / "README.md"))
    if missing:
        print("MISSING (not copied):")
        for m in missing:
            print("  ", m)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

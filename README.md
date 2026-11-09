# Symmetry-Selected Nonlinear Thermal Pumping in a Periodically Heated Mixed-Convection Channel

Repository for study 3 of the three-part series on streamwise-periodically heated
mixed convection. It contains the complete reproducibility package of the paper:
LaTeX source of both documents, the compiled PDFs, every data file produced by the
verification campaigns, the scripts that produced them, the figures, and the
quality gates that check the numbers quoted in the text.

Manuscript status: submitted to *Numerical Heat Transfer, Part B: Fundamentals*
(Taylor & Francis). The source and PDF are deposited here with the author's consent
to make every quoted figure reconstructible from the raw data.

## Contents

The paper is published as two documents: a main manuscript and an independently
numbered supplementary document whose floats, sections and equations carry the
prefix **S** (Fig. S1, Table S1, Sec. S1, Eq. (S1)). The main text runs in nine
sections -- introduction, formulation and base state, the symmetry-selected
pumping law, the exact pumping identities and their Reynolds-number dependence,
the phase-shifted two-coefficient law, numerical methods and verification,
results with the independent nonlinear DNS validation, physical interpretation
and discussion, conclusions -- and keeps the parity theorem, the two-coefficient
phase law, the reversal cartography and every key measurement. The
supplement follows the order in which a reviewer checks such a claim:
`S1` perturbative parameter maps and scaling checks, `S2` nonlinear DNS
verification and convergence, `S3` phase decomposition and term-wise balance,
`S4` physical structure of the symmetry-selected pumping mechanism, `S5` the
seven derivations and the independent shooting check, `S6` construction of the two
independent mixed kernels. Each of the 21 floats of that document is cited and
discussed from the main text, and cross-references resolve through the `xr`
package, so the two files cannot silently disagree about a number.

| Path | What it holds |
|---|---|
| `main_v2.tex` | the main manuscript (41 pages when compiled) |
| `supplementary_text.tex` | the supplementary document, S-numbered (21 pages when compiled) |
| `manuscript/main_v2.pdf` | the compiled main manuscript, build of this repository |
| `manuscript/supplementary_text.pdf` | the compiled supplementary document |
| `references.bib` | bibliography; the main manuscript carries the only reference list |
| `cover_letter_NHTB_v2.tex` | cover letter addressed to the editor (1 page) |
| `supplementary_material/` | the 64-file supplementary data package, groups A-D, with `MANIFEST.csv` (size, line count and SHA-256 of every entry) and its own `README.md` |
| `campaign_v2_symmetry_selection/` | the campaigns: raw CSV/JSON results, fitted parameters, publication figures, drivers, analysis scripts and quality gates |
| `dns_physical_validation_suite/` | the fully nonlinear DNS records: per-case measurements, the four convergence and scaling scans, the wall-flux parity test and the run report the consistency script re-derives |
| `validation_outputs_ULTIMATE_DENSE_CPU/` | the dense finite-domain critical-Peclet survey behind the endpoint classification |
| `*.png` at the root | the first-order, DNS-convergence and mechanism figures included by the two documents |
| `validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py`, `dns_physical_validation_suite.py`, `ultimate_dense_cpu_validation_publication_CORRECTED.py` at the root | the three solvers: the second-order boundary-value machinery, the fully nonlinear DNS suite, and the dense critical-Peclet survey |

## The campaigns in `campaign_v2_symmetry_selection/`

| Prefix | Campaign |
|---|---|
| `D1_*` | reversal cartography: the critical Reynolds number as a function of wavenumber, primary window and low-`k` extension |
| `D2_*` | Grashof invariance of the reversal point; `DNS_D2_grpair_*` is the nonlinear Grashof-pair run |
| `D5_*`, `C2_*`, `P2_*` | two-coefficient phase law, equal-amplitude `sin(phi)` response, direct phase sweep |
| `DNS_D345_*` | direct numerical simulation matrix: amplitude, spatial and time-step convergence, Grashof and Reynolds scans, wall-flux parity test |
| `ultimate_*` | dense finite-domain critical-Peclet survey (the data behind the endpoint classification) |
| `run_*.py`, `analyze_*.py`, `make_fig_*.py` | drivers, post-processing and figure generation |
| `restructure_sections.py`, `reorg_*.py`, `move_*.py`, `add_end_matter.py` | the section recomposition itself: each moves whole LaTeX blocks by delimiter, never by line number, so a section cannot be split mid-paragraph |
| `make_supplementary_material.py`, `make_submission_folder.py` | assemble the 64-file data package and the upload folder, with the README generated from the build logs; the folder rewrites only the two `xr` targets, to the names it ships |
| `apply_review_round2.py`, `apply_review_round2_style.py` | the referee rounds: asserted replacements of whole sentences, matched exactly or whitespace-insensitively when a paragraph was reflowed, aborting on an ambiguous match |
| `check_*.py`, `verify_manuscript_numbers.py` | the six quality gates described below |

## Reproducing

Environment: Python 3 with NumPy, SciPy and matplotlib; MiKTeX or TeX Live for the
documents. From the repository root:

```bash
py -X utf8 campaign_v2_symmetry_selection/run_campaigns.py          # perturbative campaigns
py -X utf8 campaign_v2_symmetry_selection/analyze_dns_D345.py       # DNS post-processing
py -X utf8 campaign_v2_symmetry_selection/make_supplementary_material.py
```

The DNS matrix is not re-run by these commands on its own: `run_dns_D3_D4_D5.py`
performs a dry run unless an explicit go-ahead flag is supplied, and otherwise
resumes from its own checkpoint.

To rebuild the two documents:

```bash
mkdir -p build_v2 build_cover
pdflatex -interaction=nonstopmode -output-directory=build_v2 main_v2.tex
bibtex   build_v2/main_v2
pdflatex -interaction=nonstopmode -output-directory=build_v2 main_v2.tex
pdflatex -interaction=nonstopmode -output-directory=build_v2 supplementary_text.tex
pdflatex -interaction=nonstopmode -output-directory=build_v2 supplementary_text.tex
pdflatex -interaction=nonstopmode -output-directory=build_v2 main_v2.tex
pdflatex -interaction=nonstopmode -output-directory=build_v2 supplementary_text.tex
pdflatex -interaction=nonstopmode -output-directory=build_cover cover_letter_NHTB_v2.tex
```

The passes alternate because the two documents reference each other with the `xr`
package: the supplement needs the auxiliary file of the main text, and the main text
needs the auxiliary file of the supplement for the S-numbered cross-references.
Building the supplement only at the end leaves the S-numbered cross-references
undefined in the main text, which `check_build.py` reports.

Compile from the repository root: both sources address their figures with
root-relative paths, and the `xr` blocks point at `build_v2/main_v2.aux` and
`build_v2/supplementary_text.aux`. The build directories are named as above
because the gates read their logs and auxiliary files from there.

## Quality gates

All six scripts exit non-zero on any defect, so neither document can drift away
from the data it quotes. Each resolves the repository root from its own location.

```bash
py -X utf8 campaign_v2_symmetry_selection/check_build.py                # build, labels, ASCII, prose hygiene
py -X utf8 campaign_v2_symmetry_selection/verify_manuscript_numbers.py  # every quoted number re-derived
py -X utf8 campaign_v2_symmetry_selection/check_notation.py             # notation and macro use
py -X utf8 campaign_v2_symmetry_selection/check_duplicate_text.py       # duplication scan against the companion sources
py -X utf8 campaign_v2_symmetry_selection/check_cover_letter.py         # cover-letter length and required statements
py -X utf8 campaign_v2_symmetry_selection/check_end_matter.py           # abstract length, supplementary pointers, declarations
```

`check_build.py` parses the logs and `.aux` files of **both** documents, treats them
as one corpus for dangling references and uncited floats, refuses any working
filename in the prose, requires a keyword line, requires every figure and table to be
cited in the text, requires each of the 21 supplementary floats to be cited from the
main text, checks each of the 24 graphics inputs against the file tree with its exact
case and each `xr` target against the source and auxiliary file it names, and enforces
a 45-page ceiling on the main manuscript. It needs the build
above to exist; `verify_manuscript_numbers.py` and `check_notation.py` run on the
data and the sources alone, the latter auditing both documents for the same notation
conventions.

Two gates are author-side by nature and cannot pass from this repository alone:
`check_duplicate_text.py` compares the manuscript **and** its supplementary document
against the sources of both companion studies, which are not deposited here; only the
first of them is disclosed in the submission text, and
`check_cover_letter.py` inspects the compiled cover letter, which must be built
locally from `cover_letter_NHTB_v2.tex` into a `build_cover/` directory.

## Relationship to the companion study

The paper continues the work whose first installment treats the symmetric modulation and
the non-monotonic integral framework; that study is under review at *Numerical Heat
Transfer, Part B*, and it is the **only** companion named in the manuscript and in the
cover letter, because it is the one the editor can verify. Other companion material that
is neither submitted nor published is deliberately left out of the submission text, and
the vanishing of the mean second-order component in the symmetric limit is derived in
this paper rather than quoted from unpublished work. The textual overlap is measured,
not asserted: `check_duplicate_text.py` scans the manuscript together with its
supplementary document against the companion sources and against the earlier drafts of
this paper, reports the shared twelve-word runs, and fails if prose is carried over
from a companion or if the accounting reappears in the manuscript or in the cover
letter: both documents state nothing about duplication, because the scan is a measure
taken by the author, not a claim addressed to the editor.

## Main results

- Exact parity theorem for the second-order mean thermal flux, verified in direct
  numerical simulation to a wall-flux difference of 3.1e-18 against individual fluxes
  of 1.145796985e-7 (relative 2.7e-11).
- Linearity of the pumping coefficient in the Grashof number, and Gr-independence of
  the reversal point, which survives nonlinearly: the simulated Grashof ratios are
  4.0003 and 4.0015 and the perturbative-discrepancy ratios are 3.996 and 4.015.
- A deviation model specified before the campaign, exceeded in 5 of 12 cases,
  reported and characterized.
- Two independent solvers agreeing on the same boundary-value problem to 9e-12:
  adaptive collocation and a fixed-mesh shooting method.

## License

- Text, figures and data: [CC BY 4.0](LICENSE.md)
- Code: MIT (`LICENSE_CODE.md`)

## Citation

See `CITATION.cff`. Cite the manuscript version accepted for publication; the
repository DOI reserved through this file points to whatever version you are reading.

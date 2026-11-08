# Symmetry-Selected Nonlinear Thermal Pumping in a Periodically Heated Mixed-Convection Channel

Repository for study 3 of the three-part series on streamwise-periodically heated
mixed convection. It contains the complete reproducibility package of the manuscript:
LaTeX source, the compiled PDF, every data file produced by the verification
campaigns, the scripts that produced them, the figures, and the quality gates that
check the numbers quoted in the text.

Manuscript status: submitted to *Numerical Heat Transfer, Part B: Fundamentals*
(Taylor & Francis). The source and PDF are deposited here with the author's consent
to make every quoted figure reconstructible from the raw data.

## Contents

| Path | What it holds |
|---|---|
| `main_v2.tex` | the manuscript (56 pages when compiled) |
| `manuscript/main_v2.pdf` | the compiled manuscript, build of this repository |
| `references.bib` | bibliography |
| `cover_letter_NHTB_v2.tex` | cover letter addressed to the editor |
| `supplementary_material/` | the 64-file supplementary package, groups S1-S4, with `MANIFEST.csv` (size, line count and SHA-256 of every entry) and its own `README.md` |
| `campaign_v2_symmetry_selection/` | the campaigns: raw CSV/JSON results, fitted parameters, publication figures, drivers, analysis scripts and quality gates |
| `*.png` at the root | the first-order, DNS-convergence and mechanism figures included by the manuscript |

## The campaigns in `campaign_v2_symmetry_selection/`

| Prefix | Campaign |
|---|---|
| `D1_*` | reversal cartography: the critical Reynolds number as a function of wavenumber, primary window and low-`k` extension |
| `D2_*` | Grashof invariance of the reversal point; `DNS_D2_grpair_*` is the nonlinear Grashof-pair run |
| `D5_*`, `C2_*`, `P2_*` | two-coefficient phase law, equal-amplitude `sin(phi)` response, direct phase sweep |
| `DNS_D345_*` | direct numerical simulation matrix: amplitude, spatial and time-step convergence, Grashof and Reynolds scans, wall-flux parity test |
| `ultimate_*` | dense finite-domain critical-Peclet survey (the data behind the endpoint classification) |
| `run_*.py`, `analyze_*.py`, `make_fig_*.py` | drivers, post-processing and figure generation |
| `check_*.py`, `verify_manuscript_numbers.py` | the five quality gates described below |

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

To rebuild the manuscript:

```bash
pdflatex -interaction=nonstopmode -output-directory=build main_v2.tex
bibtex   build/main_v2
pdflatex -interaction=nonstopmode -output-directory=build main_v2.tex
pdflatex -interaction=nonstopmode -output-directory=build main_v2.tex
```

Compile from the repository root: the manuscript addresses its figures with
root-relative paths.

## Quality gates

All five scripts exit non-zero on any defect, so the manuscript cannot drift away
from the data it quotes.

```bash
py -X utf8 campaign_v2_symmetry_selection/check_build.py                # build, labels, ASCII, prose hygiene
py -X utf8 campaign_v2_symmetry_selection/verify_manuscript_numbers.py  # every quoted number re-derived
py -X utf8 campaign_v2_symmetry_selection/check_notation.py             # notation and macro use
py -X utf8 campaign_v2_symmetry_selection/check_duplicate_text.py       # 12-word shingle overlap with the companion studies
py -X utf8 campaign_v2_symmetry_selection/check_cover_letter.py         # cover-letter length and required statements
```

`check_build.py` additionally refuses any working filename in the prose, requires a
keyword line, and requires every figure and table to be cited in the text.

## Relationship to the other two studies

This is one paper of a three-study series. Study 1 treats the symmetric modulation
and the non-monotonic integral framework; study 2 corrects the modal decomposition
and establishes the vanishing of the mean second-order Nusselt component. Neither is
citable yet, so both are described in the text rather than referenced, and the
present study is separated from them by the asymmetric configuration and by the
exact parity selection rule. The overlap is measured, not asserted: `check_duplicate_text.py`
compares the manuscript against the companion sources and against the earlier drafts
of this paper.

## Headline results

- Exact parity theorem for the second-order mean thermal flux, verified in direct
  numerical simulation to a wall-flux difference of 3.1e-18 against individual fluxes
  of 1.145796985e-7 (relative 2.7e-11).
- Linearity of the pumping coefficient in the Grashof number, and Gr-independence of
  the reversal point, which survives nonlinearly: the simulated Grashof ratios are
  4.0003 and 4.0015 and the perturbative-discrepancy ratios are 3.996 and 4.015.
- A pre-registered deviation model that the campaign falsified in 5 of 12 cases,
  reported as a finding rather than hidden.
- Two independent solvers agreeing on the same boundary-value problem to 9e-12:
  adaptive collocation and a fixed-mesh shooting method.

## License

- Text, figures and data: [CC BY 4.0](LICENSE.md)
- Code: MIT (`LICENSE_CODE.md`)

## Citation

See `CITATION.cff`. Cite the manuscript version accepted for publication; the
repository DOI reserved through this file points to whatever version you are reading.

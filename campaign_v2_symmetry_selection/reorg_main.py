"""Recompose the main manuscript into the nine sections requested by the review.

The moves are mechanical and block-based, never line-number based:
  * the three theory sections keep their content and take the requested titles;
  * a new Section 6, Numerical Methods and Verification, is formed by moving the
    direct-simulation formulation block (configuration, estimator, independence)
    out of the nonlinear-validation section and next to the boundary-value solver
    description;
  * the former Sections 6 (results), 7 (DNS) and 8 (mapping) plus the discussion
    subsection that reports the thirteen-case campaign merge into Section 7,
    Results and Independent Nonlinear DNS Validation;
  * the two blocks that report the campaign and its error model leave the
    discussion and join the results they describe.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "main_v2.tex"

TITLES = {
    "\\section{The Symmetry Law: Odd Parity of the\\\\ Cross-Interaction}":
        "\\section{Symmetry-Selected Nonlinear Thermal Pumping}",
    "\\section{Exact Identities and Thermal Pumping Coefficient}":
        "\\section{Exact Pumping Identities and Reynolds-Number Dependence}",
    "\\section{Exact Phase Selection and the Quadrature Pumping Coefficient}":
        "\\section{Phase-Shifted Forcing and the Two-Coefficient Law}",
    "\\section{Numerical Results and Validation}":
        "\\section{Numerical Methods and Verification}",
}


def find(lines, needle, start=0):
    for i in range(start, len(lines)):
        if needle in lines[i]:
            return i
    raise StopIteration


def block_to_next(lines, start, stops):
    """[start, end) from *start* up to the next line starting any stop."""
    for i in range(start + 1, len(lines)):
        if any(lines[i].startswith(s) for s in stops):
            return start, i
    return start, len(lines)


def main():
    lines = SRC.read_text(encoding="utf-8").splitlines(keepends=True)

    # 1. requested titles for the theory sections and for the new methods section
    for old, new in TITLES.items():
        i = find(lines, old)
        lines[i] = new + "\n"

    # 2. the DNS formulation block becomes a subsection of the methods section
    head = find(lines, "\\section{Independent Nonlinear DNS Validation}")
    span = block_to_next(lines, head, ["\\subsection{"])
    dns_block = lines[span[0]:span[1]]
    label = [t for t in dns_block if "\\label{sec:dns_validation}" in t]
    dns_block = [t for t in dns_block if t not in label]
    for j, t in enumerate(dns_block):
        if t.startswith("\\section{Independent Nonlinear DNS"):
            dns_block[j] = ("\\subsection{Direct numerical simulation: formulation, "
                            "estimator and independence}\n"
                            + label[0])
    del lines[span[0]:span[1]]

    # the results label follows the new results section, not the methods one
    m = find(lines, "\\label{sec:results}")
    lines[m] = "\\label{sec:methods}\n"
    anchor = find(lines, "\\subsection{Grashof-number scaling}")
    lines[anchor:anchor] = (dns_block
                            + ["\n", "\\clearpage\n",
                               "%---" + "-" * 62 + "\n",
                               "\\section{Results and Independent Nonlinear DNS "
                               "Validation}\n",
                               "\\label{sec:results}\n",
                               "%---" + "-" * 62 + "\n",
                               "\n"])

    # 3. the mapping section becomes the closing part of the results section
    i = find(lines, "\\section{Systematic Mapping of the Reynolds-Number Reversal}")
    lines[i] = "\\subsection{Systematic mapping of the Reynolds-number reversal}\n"

    # 4. campaign blocks move from the discussion into the results they describe
    start = find(lines, "\\subsection{Fully nonlinear DNS confirmation}")
    stop = find(lines, "\\subsection{What drives the growth of the remainder")
    campaign = lines[start:stop]
    del lines[start:stop]
    target = find(lines, "\\section{Physical Interpretation and Discussion}")
    lines[target:target] = campaign + ["%---" + "-" * 62 + "\n"]

    SRC.write_text("".join(lines), encoding="utf-8")
    for j, t in enumerate(lines, 1):
        if t.startswith(("\\section{", "\\subsection{")):
            print("%4d  %s" % (j, t.strip()[:78]))


if __name__ == "__main__":
    main()

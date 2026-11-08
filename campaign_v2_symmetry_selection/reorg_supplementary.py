"""Recompose the supplementary document into the six thematic sections S1-S6.

Moves whole LaTeX blocks by their own delimiters, never by line numbers:
  * fig:logic_diagram and tab:paper_comparison leave the phase section and open
    the physical-structure section;
  * the former physical-realization section is merged into that section, so its
    five subsections become S4.1-S4.5;
  * the derivations and the kernel construction move to the end, becoming S5
    and S6, which is the only order in which the section numbers match the
    numbering the main manuscript points at.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "supplementary_text.tex"

TITLES = {
    "Verification of the nonlinear DNS":
        "Nonlinear DNS verification and convergence",
    "Phase decomposition, term balance and transport comparison":
        "Phase decomposition and term-wise balance",
    "Two-dimensional physical structure of the mechanism":
        "Physical structure of the symmetry-selected pumping mechanism",
    "Derivations and independent checks":
        "Derivations and independent numerical checks",
    "Physical realization of the symmetry-selected pumping mechanism": None,
}


def env_span(lines, label):
    """Index span [start, end) of the float environment that owns *label*."""
    kind = "figure" if label.startswith("fig:") else "table"
    hit = next(i for i, t in enumerate(lines) if ("\\label{%s}" % label) in t)
    start = max(i for i, t in enumerate(lines[:hit])
                if t.startswith("\\begin{%s}" % kind))
    stop = min(i for i, t in enumerate(lines[start + 1:], start + 1)
               if t.startswith("\\end{%s}" % kind))
    return start, stop + 1


def cut(lines, start, stop):
    block = lines[start:stop]
    del lines[start:stop]
    return block


def find(lines, needle, start=0):
    return next(i for i, t in enumerate(lines[start:], start) if needle in t)


def main():
    lines = SRC.read_text(encoding="utf-8").splitlines(keepends=True)

    # 1. floats that belong to the physical-structure section
    logic = cut(lines, *env_span(lines, "fig:logic_diagram"))
    comp = cut(lines, *env_span(lines, "tab:paper_comparison"))

    # 2. derivations and kernel construction, which must come last
    d_start = find(lines, "\\section{Derivations")
    if lines[d_start - 1].strip() == "\\clearpage":
        d_start -= 1
    p_start = find(lines, "% ---- Physical realization")
    derivations = cut(lines, d_start, p_start)

    # 3. drop the merged header, keep its label for the new section title
    hdr = find(lines, "\\section{Physical realization")
    for i in range(hdr, hdr - 2, -1):
        if lines[i].lstrip().startswith("%"):
            hdr = i
    label_line = cut(lines, hdr, find(lines, "%====", hdr) + 1)
    label = [t for t in label_line if "\\label{" in t]

    # 4. opening prose of the merged section, moved above its field figures
    first_sub = find(lines, "\\subsection{Two-dimensional thermal structure")
    head = find(lines, "\\section{Two-dimensional physical structure")
    prose_start = find(lines, "The exact analysis above identifies")
    prose = cut(lines, prose_start, first_sub)
    while prose and not prose[0].strip():
        prose.pop(0)
    while prose and not prose[-1].strip():
        prose.pop()
    prose.append("\n")

    lines[head + 1:head + 1] = label + ["\n"] + prose + logic + ["\n"] + comp + ["\n"]

    # 5. the two derivation sections close the document
    tail = find(lines, "\\end{document}")
    lines[tail:tail] = derivations + ["\n"]

    # 6. section titles, once every block has reached its final position
    for old, new in TITLES.items():
        if new is None:
            continue
        i = find(lines, "\\section{%s}" % old)
        lines[i] = lines[i].replace(old, new)

    SRC.write_text("".join(lines), encoding="utf-8")
    for i, t in enumerate(lines, 1):
        if t.startswith(("\\section{", "\\subsection{")):
            print("%4d  %s" % (i, t.strip()[:74]))


if __name__ == "__main__":
    main()

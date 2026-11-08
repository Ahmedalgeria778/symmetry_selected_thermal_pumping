"""Add the Supplementary Material section, the declarations and the pointers.

The editorial rules need three things the manuscript did not have: a section
titled exactly "Supplementary Material" after the Conclusions that describes the
separate document, the publisher's declaration statements, and body pointers
that use the words "supplementary material" rather than an abbreviated label.
Section-level labels are added to the supplementary document so that the
description can cite its six sections by reference instead of by typed number.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "main_v2.tex"
SUPP = ROOT / "supplementary_text.tex"

SUPP_LABELS = {
    "\\section{Perturbative parameter maps": "sec:maps",
    "\\section{Nonlinear DNS verification": "sec:dns_verification",
    "\\section{Phase decomposition and term-wise": "sec:phase_balance",
    "\\section{Derivations and independent numerical": "sec:derivations",
}

END_MATTER = r"""\section*{Supplementary Material}
\begin{sloppypar}
The derivations, the verification records, the parameter maps and the field
structures that support this article are collected in a separate Supplementary
Material document, compiled independently and numbered throughout with the
prefix \textbf{S}, so that its figures, tables, equations and section numbers
cannot be confused with those of the manuscript. It contains six sections:
Sec.~\ref{sec:maps}, perturbative parameter maps and scaling checks;
Sec.~\ref{sec:dns_verification}, nonlinear DNS verification and convergence,
which carries the amplitude, grid and time-step sequences
Sec.~\ref{subsec:dns_sequences} summarises; Sec.~\ref{sec:phase_balance}, phase
decomposition and term-wise balance; Sec.~\ref{sec:physical_mechanism}, physical
structure of the symmetry-selected pumping mechanism, which includes the spatial
rectification analysis of Sec.~\ref{sec:spatial_rectification};
Sec.~\ref{sec:derivations}, derivations and independent numerical checks; and
Sec.~\ref{sec:kernels}, construction of the two independent mixed kernels. Every
figure and table of that document is cited from the manuscript, and the pointers
in the text read "Sec.~S\#.\# of the supplementary material".
\end{sloppypar}

"""

DECLARATIONS = r"""\section*{Disclosure statement}
The author reports no potential conflict of interest.

\section*{Author contributions}
The author conceived and designed the study, derived the theory, implemented the
numerical methods, ran every verification and direct-numerical-simulation
campaign, processed the data, prepared the figures, tables and the data package,
and wrote and edited the manuscript.

\section*{Funding}
This research did not receive any specific grant from funding agencies in the
public, commercial, or not-for-profit sectors.

"""


def main():
    supp = SUPP.read_text(encoding="utf-8")
    for prefix, label in SUPP_LABELS.items():
        if ("\\label{%s}" % label) in supp:
            continue
        supp = re.sub("(?m)^%s[^\n]*\n" % re.escape(prefix),
                      lambda m: m.group(0) + "\\label{%s}\n" % label,
                      supp, count=1)
    SUPP.write_text(supp, encoding="utf-8")

    text = MAIN.read_text(encoding="utf-8")
    text = re.sub(r"Supplementary Sec\.~\\ref\{([^}]*)\}",
                  r"Sec.~\\ref{\1} of the supplementary material", text)

    anchor = "\\begin{acknowledgments}"
    text = text.replace(anchor, END_MATTER + anchor, 1)
    anchor = "\\section*{Data availability}"
    text = text.replace(anchor, DECLARATIONS + anchor, 1)
    MAIN.write_text(text, encoding="utf-8")

    print("pointers reworded: %d" %
          len(re.findall(r"of the supplementary material", text)))
    for title in ("Supplementary Material", "Disclosure statement",
                  "Author contributions", "Funding", "Data availability"):
        print("  \\section*{%s}: %s" % (title,
                                        "\\section*{%s}" % title in text))


if __name__ == "__main__":
    main()

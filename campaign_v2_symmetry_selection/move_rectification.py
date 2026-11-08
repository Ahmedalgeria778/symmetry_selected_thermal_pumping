"""Move the detailed spatial-rectification analysis from the discussion to S4.

Review move C asks the main manuscript to keep a mechanism explained in a few
paragraphs and to send the detailed spatial analysis to the supplementary
document. The moved prose is transferred verbatim: nothing is re-derived here,
and the numerical statements stay exactly where their data are.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "main_v2.tex"
SUPP = ROOT / "supplementary_text.tex"

CONDENSED = r"""\subsection{Mechanism: buoyancy-mediated spatial rectification}
\label{sec:spatial_rectification}

The mechanism is a two-stage process. The boundary asymmetry creates simultaneous
symmetric and antisymmetric thermal modes; buoyancy converts them into transverse
velocity modes, and the resulting velocity--temperature products generate the odd
second-order forcing, whose magnitude carries the linear factor $\Grr$. The content
of that forcing is most compactly stated as a pair of moment conditions,
\begin{equation}
  \int_{-1}^{1}\mathcal{E}_{as}(Y)\,dY=0,
  \qquad
  \int_{-1}^{1}Y\,\mathcal{E}_{as}(Y)\,dY\neq0,
  \label{eq:moment_pair}
\end{equation}
which is what makes the channel a \emph{transverse spatial rectifier}: the zeroth
moment cancels by parity, so no net mean thermal energy is generated, while the first
transverse moment survives and is converted directly into the measurable wall-flux
difference by the Green identity Eq.~\eqref{eq:PiT_integral}. The selection theorem
fixes which wall-flux combinations are admissible; the coupled transport problem fixes
how strongly and in which direction the admissible redistribution occurs. That
division of labour is why the reversal near $\Rey_c\simeq30.84$ at fixed $\Pe$
reports a change of the hydrodynamic pathway and not of the thermal forcing class.
Two limits carry the same reading: the $O(k^{2})$ suppression of the rectified moment
as the wall pattern becomes slowly varying, and the finite-domain endpoint of the
zero-crossing branch, where the zero set leaves the admissible domain rather than
losing stability. The field maps, the secondary circulation and the detailed spatial
analysis of both limits are developed in
Supplementary Sec.~\ref{subsec:rectification}.
"""


def find(lines, needle, start=0):
    for i in range(start, len(lines)):
        if needle in lines[i]:
            return i
    raise StopIteration


def main():
    lines = MAIN.read_text(encoding="utf-8").splitlines(keepends=True)

    start = find(lines, "\\subsection{Buoyancy as the hydrodynamic mediator}")
    stop = find(lines, "\\subsection{What is exact and what is dynamical?}")
    moved = lines[start:stop]
    lines[start:stop] = [CONDENSED + "\n"]
    MAIN.write_text("".join(lines), encoding="utf-8")

    # the parity-realization subsection that precedes it belongs to the same text
    keep = [t for t in moved if t.startswith("\\subsection")]

    body = []
    for t in moved:
        if t.startswith("\\subsection{Physical realization"):
            continue
        if t.startswith("\\subsection{Buoyancy as the hydrodynamic mediator}"):
            body.append("\\subsection{Spatial rectification and parameter-space "
                        "structure}\n\\label{subsec:rectification}\n\n")
            body.append("The buoyancy-mediated reading of the previous sections is a "
                        "two-stage process: the boundary asymmetry creates simultaneous "
                        "symmetric and antisymmetric thermal modes, buoyancy converts them "
                        "into transverse velocity modes, and the resulting "
                        "velocity--temperature products generate the odd second-order "
                        "forcing whose magnitude carries the linear factor $\\Grr$. "
                        "Table~\\ref{tab:paper_comparison} states the contrast with the "
                        "symmetry-preserving forcing term by term.\n\n")
            continue
        if t.startswith("\\subsection{Spatial rectification"):
            continue
        if "\\label{sec:spatial_rectification}" in t:
            continue
        if t.startswith("%===="):
            continue
        body.append(t)
    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    body.append("\n")

    supp = SUPP.read_text(encoding="utf-8").splitlines(keepends=True)
    cut = find(supp, "\\clearpage")
    supp[cut:cut] = ["\n"] + body
    SUPP.write_text("".join(supp), encoding="utf-8")
    print("moved %d lines into Supplementary Sec.~S4; kept %d condensed lines" %
          (len(body), CONDENSED.count("\n")))
    print("subsections carried by the moved block:", keep)


if __name__ == "__main__":
    main()

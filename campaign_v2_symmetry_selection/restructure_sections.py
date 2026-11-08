"""Plan item E: collapse the manuscript to the recommended 11-section structure.

The two trailing sections -- "Physical realization of the symmetry-selected
pumping mechanism" and "Discussion: spatial rectification and parameter-space
structure" -- become subsections of one interpretation-and-discussion section,
and the closing subsections (exact vs dynamical, limitations, generalization,
outlook) move to the end so that the outlook is the last thing before the
conclusions.

Every split is asserted against an exact heading string, so the script fails
loudly on a manuscript that has drifted. A one-time backup is written next to
the source.
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, "main_v2.tex")
BAK = os.path.join(ROOT, "main_v2.tex.bak-preE")

DISC_HEAD = "\\section{Discussion and Physical Interpretation}"
DISC_NEW = "\\section{Physical Interpretation and Discussion}\n\\label{sec:discussion}"
PHYS_HEAD = "\\section{Physical realization of the symmetry-selected pumping mechanism}"
PHYS_NEW = ("\\subsection{Physical realization of the symmetry-selected pumping "
            "mechanism}")
SPAT_HEAD = "\\section{Discussion: spatial rectification and parameter-space structure}"
SPAT_NEW = "\\subsection{Spatial rectification and parameter-space structure}"
CONC_HEAD = "\\section{Conclusions}"

PHASE_SUB = "\\subsection{Quantitative phase analysis of the reversal}"
MED_SUB = "\\subsection{Buoyancy as the hydrodynamic mediator}"
EXACT_SUB = "\\subsection{What is exact and what is dynamical?}"


def once(lines, needle, what):
    hits = [i for i, l in enumerate(lines) if l.strip() == needle]
    if len(hits) != 1:
        sys.exit("ABORT: expected exactly one %s, found %d" % (what, len(hits)))
    return hits[0]


def trim(block):
    """Drop decorative comment rules and blank lines at both ends."""
    while block and (not block[0].strip() or block[0].lstrip().startswith("%")):
        block = block[1:]
    while block and (not block[-1].strip() or block[-1].lstrip().startswith("%")):
        block = block[:-1]
    return block


def main():
    text = io.open(TEX, encoding="utf-8").read()
    lines = text.split("\n")

    i_disc = once(lines, DISC_HEAD, "discussion section heading")
    i_phys = once(lines, PHYS_HEAD, "physical-realization section heading")
    i_spat = once(lines, SPAT_HEAD, "spatial-rectification section heading")
    i_conc = once(lines, CONC_HEAD, "conclusions heading")
    i_phase = once(lines, PHASE_SUB, "phase-analysis subsection")
    i_med = once(lines, MED_SUB, "buoyancy-mediator subsection")
    i_exact = once(lines, EXACT_SUB, "exact-vs-dynamical subsection")
    if not (i_disc < i_phase < i_med < i_exact < i_phys < i_spat < i_conc):
        sys.exit("ABORT: headings are not in the expected order")

    if not os.path.exists(BAK):
        io.open(BAK, "w", encoding="utf-8", newline="\n").write(text)
        print("backup written: %s" % BAK)

    pre = lines[:i_disc - 1]                                  # keep the rule style
    intro = trim(lines[i_disc + 2:i_phase])                  # section opening prose
    phase = lines[i_phase:i_med]
    mediator = lines[i_med:i_exact]
    closing = trim(lines[i_exact:i_phys])
    phys = trim(lines[i_phys:i_spat])
    spatial = trim(lines[i_spat:i_conc])
    tail = lines[i_conc:]

    if not phys[0].strip().startswith("\\section"):
        sys.exit("ABORT: physical block does not start on its heading")
    phys[0] = PHYS_NEW
    n_demoted = 0
    for i in range(1, len(phys)):
        if phys[i].startswith("\\subsection{"):
            phys[i] = "\\subsubsection{" + phys[i][len("\\subsection{"):]
            n_demoted += 1
    spatial[0] = SPAT_NEW

    rule = "%--------------------------------------------------------------"
    out = (pre
           + [lines[i_disc - 1], DISC_NEW, rule, ""]
           + intro + [""]
           + phase + phys + mediator + spatial + closing
           + [""] + tail)

    io.open(TEX, "w", encoding="utf-8", newline="\n").write("\n".join(out))

    print("physical-realization subsections demoted: %d" % n_demoted)
    sections = re.findall(r"^\\section\{(.*)\}", "\n".join(out), re.M)
    print("\nsections now (%d):" % len(sections))
    for i, s in enumerate(sections, 1):
        print("  %2d. %s" % (i, s))
    return 0


if __name__ == "__main__":
    sys.exit(main())

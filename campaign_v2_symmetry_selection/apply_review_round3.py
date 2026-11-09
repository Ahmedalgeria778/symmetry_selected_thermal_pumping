#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Apply the third referee round: the abstract, the manuscript facts, the letter, the supplement.

Every replacement is asserted against the source, with a whitespace-insensitive
retry so a reflowed paragraph cannot skip a fix. The abstract is installed as a
block because it is replaced whole.

    py -X utf8 campaign_v2_symmetry_selection/apply_review_round3.py
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(name):
    return io.open(os.path.join(ROOT, name), encoding="utf-8").read()


def apply(text, old, new, path, count=1):
    n = text.count(old)
    if n == 0:
        pat = re.compile(r"\s+".join(re.escape(w) for w in old.split()))
        hits = list(pat.finditer(text))
        if not hits:
            raise SystemExit("NOT FOUND in %s:\n%s" % (path, old[:220]))
        if len(hits) != count:
            raise SystemExit("AMBIGUOUS (%d) in %s:\n%s"
                             % (len(hits), path, old[:220]))
        print("  reflowed match in %s: %s..." % (path, old[:60].replace("\n", " ")))
        start, stop = hits[0].span()
        return text[:start] + new + text[stop:]
    if n != count:
        raise SystemExit("AMBIGUOUS (%d) in %s:\n%s" % (n, path, old[:220]))
    return text.replace(old, new)


ABSTRACT = r"""\begin{abstract}
Nonlinear redistribution of the wall heat flux in a laminar mixed-convection
channel with spatially periodic wall heating of unequal amplitude obeys an exact
symmetry-selection rule. Decomposing the first-order thermal response into
symmetric and antisymmetric components, we prove their mixed quadratic
interaction is odd across the channel, so second-order mean wall-flux
corrections are equal and opposite; symmetric self-interactions produce no
directional redistribution at that order. The parity law defines a
thermal-pumping coefficient obeying an exact Green-function identity, from which
follow its linearity in the Grashof number, its quadratic suppression at long
wavelength, and a zero independent of the buoyancy amplitude. At fixed P\'eclet
number the Reynolds number controls the pumping direction: the reversal locus,
mapped at seven wavenumbers, shifts from $\Rey=79.1$ at $k=0.4$ to $\Rey=8.5$ at
$k=2.0$, and a term-wise decomposition shows a magnitude race between
streamwise-advection and wall-normal terms, not a phase flip.
Phase-shifted forcing adds an independent quadrature coefficient that can
reverse the pumping alone; at equal amplitudes the response is
exactly proportional to $\sin\phi$. A fully nonlinear direct numerical
simulation, independent of the perturbative solver, matches that coefficient to
$1.52\%$, returns the predicted reversal sign at four wavenumbers,
shows parity cancellation at $2.7\times10^{-11}$ relative and Grashof scaling at
finite amplitude. In five of twelve non-degenerate cases the deviation exceeds
the bound specified before the campaign, independent of amplitude and growing
with wavenumber. Symmetry selects the structure of the nonlinear
redistribution; buoyancy-mediated dynamics select its amplitude, its sign, and
the control parameter that governs it.
\end{abstract}"""


def replace_abstract(text):
    pat = re.compile(r"\\begin\{abstract\}.*?\\end\{abstract\}", re.S)
    hits = list(pat.finditer(text))
    if len(hits) != 1:
        raise SystemExit("abstract block found %d times" % len(hits))
    start, stop = hits[0].span()
    return text[:start] + ABSTRACT + text[stop:]


MAIN = "main_v2.tex"
SUPP = "supplementary_text.tex"
LETTER = "cover_letter_NHTB_v2.tex"

edits_main = [
    # 1. Table 5 caption: the table has no D3b row (the k=0.6 pair is labelled
    #    D3 there) and the ladder run carries its own conversion factor.
    (r"""for the D3, D3b and D4 families,""",
     r"""for the D3 and D4 families,""", 1),
    (r"""quadrature runs; the ladder flux gap of""",
     r"""quadrature runs; the ladder deviation of""", 1),

    # 4. five of the twelve runs exceed the bound, not two.
    (r"""two places where its deviation exceeds the bound specified before the campaign.""",
     r"""five cases where its deviation exceeds the bound specified before the campaign.""", 1),

    # 7. "cost" is lab shorthand.
    (r"""cost $491$\,s in total""", r"""took $491$\,s in total""", 1),

    # 5. the summary table records levels only; the mesh, time step and
    #    integration length are recorded in the supplement.
    (r"""Table~\ref{tab:dns_validation_summary} ($N_x=64$""",
     r"""Table~\ref{tab:dns_validation_summary} and in Sec.~\ref{subsec:dns_sequences} of the
supplementary material ($N_x=64$""", 1),

    # 2. Re=20 and Re=30 sit on the same side of Re_c=30.838, and 1.29e-7 is a
    #    deviation (predicted minus measured), not a flux level.
    (r"""Raising $\varepsilon_a$ fourfold, to $0.02$, moves the flux gap from the
  $1.29\times10^{-7}$ measured at the reference amplitude on the other side of
  the same crossing ($k=0.8$, $\Rey=20$) to $3.25\times10^{-7}$""",
     r"""Raising $\varepsilon_a$ fourfold, to $0.02$, moves the deviation from the
  $1.29\times10^{-7}$ measured at the reference amplitude at the same wavenumber
  ($k=0.8$, $\Rey=20$) to $3.25\times10^{-7}$""", 1),

    # 3. in coefficient space the deviation falls to the floor, so "does not
    #    improve" contradicted the numbers two sentences above.
    (r"""The deviation is
  therefore \emph{not} absolute in $q_+-q_-$ and does not improve with the
  forcing amplitude: the coefficient-space reading is the correct one, and a""",
     r"""The deviation is
  therefore \emph{not} an absolute constant in $q_+-q_-$: it grows with the
  forcing amplitude in flux space while remaining pinned at the coefficient-space
  floor, so the coefficient-space reading is the correct one, and a""", 1),

    # 6. 1.42e-4 against a floor of 2.71e-4 is about half, not half.
    (r"""deviation, $1.42\times10^{-4}$, is half the pre-specified floor""",
     r"""deviation, $1.42\times10^{-4}$, is about half the pre-specified floor""", 1),
]

edits_supp = [
    # 1. the four-controls sentence ran on with a dash and a comma splice.
    (r"""figures tighten one numerical control at a time around the reference case---the
wall amplitude, the streamwise and transverse resolution, the time step and the
Grashof pair, and the fifth, with the table that closes the section, sweeps the
Reynolds number across the reversal.""",
     r"""figures tighten one numerical control at a time around the reference case: the
wall amplitude, the streamwise and transverse resolution, the time step, and the
Grashof pair. The fifth, with the table that closes the section, sweeps the
Reynolds number across the reversal.""", 1),

    # 2. the sequence is not monotone because two levels round to the same
    #    tabulated value; say that instead of calling the change "extremely small".
    (r"""The extremely small change in $\PiT$ over this sequence shows that the reported""",
     r"""The differences at $2.0\times10^{-3}$ and $1.0\times10^{-3}$ are at the level of
round-off of the tabulated values, and the sequence is stationary to better than
$2\times10^{-6}$ throughout. That change in $\PiT$ over the sequence shows that the reported""", 1),
]

edits_letter = [
    # 1. at equal amplitudes the in-phase coefficient vanishes; only the
    #    quadrature term survives, so there is no sum of two coefficients.
    (r"""Phase-shifted wall modulation gives the response two pumping coefficients
of opposite sign whose sum is proportional to $\sin\phi$ at equal amplitudes.""",
     r"""Phase-shifted wall modulation adds a quadrature coefficient; at equal wall
amplitudes the response is proportional to $\sin\phi$, with no higher harmonic.""", 1),

    # 2. the manuscript says "reproduced"; one word for the same claim.
    (r"""or derivation is shared between the two papers.""",
     r"""or derivation of that study is reproduced here.""", 1),

    # 3 and 4.
    (r"""campaign scripts, raw tables and figure generators are deposited publicly at""",
     r"""campaign scripts, raw tables and figure-generation scripts are deposited publicly at""", 1),
    (r"""symmetry_selected_thermal_pumping}; please judge that
deposit against the journal's preprint policy.""",
     r"""symmetry_selected_thermal_pumping}. The deposit is
available for the editor's reference.""", 1),
]


def main():
    if "--abstract" in sys.argv[1:]:
        print("--- %s: abstract only" % MAIN)
        text = replace_abstract(load(MAIN))
        io.open(os.path.join(ROOT, MAIN), "w", encoding="utf-8",
                newline="\n").write(text)
        print("  written")
        return 0

    print("--- %s" % MAIN)
    text = replace_abstract(load(MAIN))
    for old, new, count in edits_main:
        text = apply(text, old, new, MAIN, count)
    io.open(os.path.join(ROOT, MAIN), "w", encoding="utf-8", newline="\n").write(text)
    print("  written")

    for name, edits in ((SUPP, edits_supp), (LETTER, edits_letter)):
        print("--- %s" % name)
        text = load(name)
        for old, new, count in edits:
            text = apply(text, old, new, name, count)
        io.open(os.path.join(ROOT, name), "w", encoding="utf-8", newline="\n").write(text)
        print("  written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

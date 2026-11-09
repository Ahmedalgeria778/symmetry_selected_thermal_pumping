#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Apply the second referee round: factual corrections and the prose de-escalation.

Every replacement is asserted. A match that fails on the exact string is retried
whitespace-insensitively, so a reflowed paragraph cannot silently skip a fix.

    py -X utf8 campaign_v2_symmetry_selection/apply_review_round2.py
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


MAIN = "main_v2.tex"
SUPP = "supplementary_text.tex"

# ----------------------------------------------------------------- main text
edits_main = [
    # item 1: one channel, everywhere
    (r"\title{\textbf{Symmetry-Selected Nonlinear Thermal Pumping in Periodically Heated Mixed-Convection Channels}}",
     r"\title{\textbf{Symmetry-Selected Nonlinear Thermal Pumping in a Periodically Heated Mixed-Convection Channel}}", 1),

    # item 2: the abstract quoted a decimal count that does not match 3.1e-18;
    # give the precision the way Sections 5 and 7 state it, and de-escalate
    # "falsifies".
    ("""A fully nonlinear direct numerical simulation converges to
within 1.52\\% of the perturbative coefficient, returns the predicted reversal
sign on both sides of the crossing at four wavenumbers, cancels the parity
partner to nineteen decimals in the wall-flux difference, and tests the
Grashof scaling at finite amplitude. The same thirteen-case campaign falsifies
part of its own pre-registered error model: five of the twelve non-degenerate
cases exceed the bound, by a deviation independent of the forcing amplitude and
growing with the wavenumber.""",
     """A fully nonlinear direct numerical simulation converges to
within 1.52\\% of the perturbative coefficient, returns the predicted reversal
sign on both sides of the crossing at four wavenumbers, cancels the parity
partner in the wall-flux difference to eleven digits relative to the flux
itself, and tests the Grashof scaling at finite amplitude. In five of the
twelve non-degenerate cases of the same thirteen-case campaign the deviation
exceeds the bound specified before the campaign, by an amount independent of
the forcing amplitude and growing with the wavenumber.""", 1),

    # item 7: equal amplitudes only forbid the pumping when they are in phase
    ("""For identical upper and lower wall modulations, $\\varepsilon_1=\\varepsilon_2$, the leading second-order
wall-flux redistribution is symmetry-forbidden.""",
     """For identical in-phase upper and lower wall modulations,
$\\varepsilon_1=\\varepsilon_2$ with $\\phi=0$, the leading second-order
wall-flux redistribution is symmetry-forbidden.""", 1),

    # item 6: nothing of this study is published yet
    ("""quadrature coefficient $\\PiT^{\\perp}$ has the \\emph{opposite sign} to the
  published in-phase coefficient $\\PiT$.""",
     """quadrature coefficient $\\PiT^{\\perp}$ has the \\emph{opposite sign} to the
  in-phase coefficient $\\PiT$ determined here.""", 1),

    # items 6 and 10: the comparison is with the direct solve of the same case,
    # and the old pointer named the section the sentence sits in
    ("""The in-phase value reproduces the authoritative
single-case calculation of Section~\\ref{sec:phase} to a relative discrepancy of
$4.3\\times10^{-13}$, which validates the kernel decomposition
\\eqref{eq:mixed_forcing_two} against the published solver.""",
     """The in-phase value agrees with the direct single-case solve of the
reference state to a relative discrepancy of $4.3\\times10^{-13}$, which
validates the kernel decomposition \\eqref{eq:mixed_forcing_two} against the
reference solver.""", 1),

    # item 9: the constraint is on the flow, not on a flux
    ("""  Zero-bulk-flux residual & $1.1\\times10^{-16}$ \\\\""",
     """  Zero-bulk-flow residual & $1.1\\times10^{-16}$ \\\\""", 1),

    # item 4: the k=0.8 crossing is a perturbative quantity, and the quoted
    # digits did not match the archived value
    ("""  sampled. The $k=0.8$ value reproduces the independent determination from the
  DNS-based sweep, $\\Rey_c=30.838496518$, to nine significant digits.}""",
     """  sampled. The $k=0.8$ entry is the value used to place the DNS cases around
  their own crossing (Sec.~\\ref{sec:outlook}); it agrees with the four
  Grashof-number determinations at that wavenumber, $\\Rey_c=30.838496520$, to
  the displayed digits, and the crossing location itself remains a perturbative
  quantity.}""", 1),

    # item 5: 3.25e-7 against 1.29e-7 is a factor 2.5, not a factor four
    ("""Raising $\\varepsilon_a$ fourfold
  multiplies the flux gap, $3.25\\times10^{-7}$ against the
  $1.29\\times10^{-7}$ measured at the reference amplitude and the nearest
  bracket point of the same wavenumber ($k=0.8$, $\\Rey=20$), but leaves the
  coefficient-space deviation at
  $2.709\\times10^{-4}$---that is $1.00$ times the absolute floor---and the
  measured sign is positive where the prediction is negative.""",
     """Raising $\\varepsilon_a$ fourfold, to $0.02$, moves the flux gap from the
  $1.29\\times10^{-7}$ measured at the reference amplitude on the other side of
  the same crossing ($k=0.8$, $\\Rey=20$) to $3.25\\times10^{-7}$: a factor of
  $2.5$, not of $4$. The coefficient-space deviation is nevertheless left at
  $2.709\\times10^{-4}$, one time the absolute floor, and the measured sign is
  positive where the prediction is negative.""", 1),

    # style: the adjective "authoritative" carried no information
    ("""Brent iteration on the authoritative perturbative solver for four Grashof numbers:""",
     """Brent iteration on the perturbative solver for four Grashof numbers:""", 1),
    ("""map the reversal locus $\\Rey_c(k)$ with the authoritative perturbative solver over""",
     """map the reversal locus $\\Rey_c(k)$ with the perturbative solver over""", 1),
    ("""from the same authoritative solver.""", """from the same solver.""", 1),
    ("""These tests support the use of the adaptive collocation result as the authoritative CPU reference while the""",
     """These tests support the use of the adaptive collocation result as the reference CPU solution while the""", 1),
    ("""accompany---but the spatial discretisation is not: on the published convergence""",
     """accompany, but the spatial discretisation is not: on the convergence""", 1),

    # style: the §7 lead paragraph described its own organisation
    ("""The evidence is ordered as the claims are. The perturbative scalings come first,
because within the model they are exact and they fix what the nonlinear runs are
supposed to reproduce. The reversal map follows, because it is the prediction the
simulation is asked to bracket. The campaign that tests the whole chain at finite
amplitude closes the section, together with the two places where it falsifies its
own pre-registered error model.""",
     """The perturbative scalings come first, because within the model they are exact
and they fix what the nonlinear runs have to reproduce. The reversal map follows,
because it is the prediction the simulation is asked to bracket. The campaign
that tests the whole chain at finite amplitude closes the section, including the
two places where its deviation exceeds the bound specified before the campaign.""", 1),

    # style: remove the superlative, keep the two measurements
    ("""  DNS at a relative precision of $2.7\\times10^{-11}$---eleven digits beyond the
  flux itself. This is the sharpest single confirmation in the campaign, and it
  is independent of the symmetric-forcing case of
  Table~\\ref{tab:dns_validation_summary} at $8\\times10^{-18}$.""",
     """  DNS at a relative precision of $2.7\\times10^{-11}$, eleven digits relative to
  the flux itself. This test is independent of the symmetric-forcing case of
  Table~\\ref{tab:dns_validation_summary}, which gives $8\\times10^{-18}$.""", 1),

    # style: neutral subsection titles, sober statement of the exceedance
    ("""\\subsection{What the campaign changes about the error model}""",
     """\\subsection{Deviation from the pre-specified bound}""", 1),
    ("""The campaign falsifies part of its own pre-registered error model. Of the twelve
non-degenerate cases, seven satisfy""",
     """The bound was exceeded in five of the twelve non-degenerate cases; we report and
characterize the excess. Seven cases satisfy""", 1),
    ("""\\subsection{What is exact and what is dynamical?}""",
     """\\subsection{Exact and dynamical statements}""", 1),
    ("""\\subsection{What drives the growth of the remainder with $k$}""",
     """\\subsection{Growth of the remainder with wavenumber}""", 1),
    ("""\\emph{What symmetry selects.} The first-order thermal fields decompose uniquely""",
     """\\emph{Symmetry-selected structure.} The first-order thermal fields decompose uniquely""", 1),
    ("""\\emph{What the dynamics select.} The redistribution is measured by the""",
     """\\emph{Dynamically selected amplitude and sign.} The redistribution is measured by the""", 1),
    ("""This is a selection rule, not a suggestion: it fixes which
combinations of wall forcing can pump and in which quadrature.""",
     """The rule fixes which combinations of wall forcing can pump, and in which
quadrature.""", 1),
    ("""The campaign also falsifies part of its own pre-registered error
model: five of the twelve non-degenerate cases exceed $\\delta_\\Pi$, with an
excess that is independent of the forcing amplitude and grows with the
wavenumber""",
     """The bound specified before the campaign is exceeded
by five of the twelve non-degenerate cases, with an excess that is independent
of the forcing amplitude and grows with the wavenumber""", 1),

    # style: boxes reserved for the parity theorem and the pumping identities
    ("""\\boxed{\\widehat U^{(j)}=-\\frac{i}{k}(\\widehat V^{(j)})'.}""",
     """\\widehat U^{(j)}=-\\frac{i}{k}(\\widehat V^{(j)})'.""", 1),
    ("""\\boxed{\\mathcal E_{as}(-Y)=-\\mathcal E_{as}(Y)}.""",
     """\\mathcal E_{as}(-Y)=-\\mathcal E_{as}(Y).""", 1),
    ("""  \\boxed{
  \\Pe_c=10.871940554611,\\qquad k_c=10,
  }""",
     """  \\Pe_c=10.871940554611,\\qquad k_c=10,""", 1),
]

# ------------------------------------------------------------ supplement
edits_supp = [
    ("""\\title{\\textbf{Supplementary Material for ``Symmetry-Selected Nonlinear Thermal
       Pumping in Periodically Heated Mixed-Convection Channels''}}""",
     """\\title{\\textbf{Supplementary Material for ``Symmetry-Selected Nonlinear Thermal
       Pumping in a Periodically Heated Mixed-Convection Channel''}}""", 1),
    ("""Grashof pair---and the fifth, with the table that closes the section, sweeps the
Reynolds number across the reversal. They fix the error levels a reader needs in
order to judge the DNS, and they are the detail from which the manuscript keeps
only the headline numbers.""",
     """Grashof pair, and the fifth, with the table that closes the section, sweeps the
Reynolds number across the reversal. They fix the error levels a reader needs in
order to judge the DNS; the manuscript quotes only the summary values.""", 1),
]

for path, edits in ((MAIN, edits_main), (SUPP, edits_supp)):
    text = load(path)
    print("--- " + path)
    for old, new, count in edits:
        text = apply(text, old, new, path, count)
    text = text.replace("pre-registered", "pre-specified")
    io.open(os.path.join(ROOT, path), "w", encoding="utf-8",
            newline="\n").write(text)
    print("  written")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Apply the second referee round, part two: the prose pass.

The reviewer's complaint is a writing register that reads as assisted output:
dual-adjective contrasts ("exact rather than rhetorical"), meta sentences that
announce the design ("the design is not accidental"), the "not X but Y" figure,
boxes around equations that are not results, and "exact" used for what is really
a measured agreement. Every replacement is asserted against the source, with a
whitespace-insensitive retry so a reflowed paragraph cannot skip a fix.

    py -X utf8 campaign_v2_symmetry_selection/apply_review_round2_style.py
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

edits_main = [
    # --- boxes only where a result is stated: keep the theorem, the Green
    #     identity and Proposition 2, drop them from the corollary, the
    #     long-wavelength limit and the restatement in the Conclusions.
    (r"""\boxed{\PiT(k,\Pe,\Grr;\Rey,\Pr)=\Grr A(k,\Rey,\Pr).}""",
     r"""\PiT(k,\Pe,\Grr;\Rey,\Pr)=\Grr A(k,\Rey,\Pr).""", 1),
    (r"""\boxed{A(k,\Rey,\Pr)=A_2(\Rey,\Pr)k^2+o(k^2),\qquad k\to0.}""",
     r"""A(k,\Rey,\Pr)=A_2(\Rey,\Pr)k^2+o(k^2),\qquad k\to0.""", 1),
    ("""\\boxed{\n\\Nu_{2,0}(+1) = -\\Nu_{2,0}(-1)\n}""",
     """\\Nu_{2,0}(+1) = -\\Nu_{2,0}(-1)""", 1),

    # --- the dual-adjective and meta tics
    ("""that the division is exact rather than rhetorical: the first statement holds as a""",
     """that the division is a difference in kind: the first statement holds as a""", 1),
    ("""which is what makes the parity assignment unique rather than accidental, is spelled out""",
     """which is what fixes the parity assignment, is spelled out""", 1),
    ("""This corollary is a structural prediction rather than an asymptotic one, and it is directly testable.""",
     """This corollary is a structural prediction, and it is directly testable.""", 1),
    ("""The DNS is consequently an independent test of the perturbative mechanism rather than a second implementation of the same boundary-value calculation.""",
     """The DNS is therefore an independent test of the perturbative mechanism: it reuses none of the boundary-value calculation.""", 1),
    ("""confirms that the transition is transversal rather than a sampling artefact: with""",
     """confirms that the transition is transversal: with""", 1),
    ("""the nonlinear evidence that $\\PiT$ is the leading transport coefficient of the symmetric--antisymmetric interaction rather than an amplitude artefact.""",
     """the nonlinear evidence that $\\PiT$ is the leading transport coefficient of the symmetric--antisymmetric interaction.""", 1),
    ("""First, the reversal is systematic rather than isolated: it""",
     """First, the reversal is systematic: it""", 1),
    ("""each wavenumber is tested on both sides of its own crossing rather than on a common Reynolds interval: at""",
     """each wavenumber is tested on both sides of its own crossing: at""", 1),
    ("""The pattern is systematic rather than random:""",
     """The pattern is systematic:""", 1),
    ("""The design is not accidental: candidate Reynolds numbers closer to the crossing""",
     """Candidate Reynolds numbers closer to the crossing""", 1),
    ("""the transition concerns the geometry of the zero set in a finite parameter domain rather than""",
     """the transition concerns the geometry of the zero set in a finite parameter domain, and not""", 1),
    ("""That test is identified as\nfuture work rather than asserted here.""",
     """That test is left for future work.""", 1),
    ("""averaging leaves two independent invariants rather than one. The wall flux""",
     """averaging leaves two independent invariants. The wall flux""", 1),

    # --- the "not X but Y" figure
    ("""interest: it is the asymmetric configuration, the exact parity selection rule that follows from""",
     """interest: it is the asymmetric configuration, the parity selection rule that follows from""", 1),
    ("""it, and the dynamical consequences of that rule. The contrast is not one of parameters but of\nwhich quantity is observable at all.""",
     """it, and the dynamical consequences of that rule. The contrast lies in which quantity is\nobservable at all, not in the parameter values.""", 1),
    ("""reading is therefore not that the two walls stop exchanging heat, but that the nonlinear\ninteraction is symmetry-forbidden from rectifying into a net transport: the pumping is\nselected, not suppressed. Buoyancy then enters only as the amplitude mediator, which is""",
     """reading is that the two walls keep exchanging heat while the nonlinear interaction is\nsymmetry-forbidden from rectifying into a net transport. Buoyancy enters only as the\namplitude mediator, which is""", 1),

    # --- em dashes used as parentheses in prose
    ("""periodic energy balance, to zero --- the endpoint representation of Supplementary""",
     """periodic energy balance, to zero; the endpoint representation of Supplementary""", 1),
    ("""fluxes --- so that in the symmetric limit the response measured here""",
     """fluxes, so that in the symmetric limit the response measured here""", 1),
    ("""of the observed zero-crossing branch for $(\\Rey,\\Grr)=(10,100)$ on $0.02\\le k\\le10$ --- not\na universal critical P\\'eclet number, not an instability threshold.""",
     """of the observed zero-crossing branch for $(\\Rey,\\Grr)=(10,100)$ on $0.02\\le k\\le10$, that is\nneither a universal critical P\\'eclet number nor an instability threshold.""", 1),
    ("""The residual departures of the two measured ratios from exactly four---$0.1\\%$\nand $0.4\\%$---are not interpreted.""",
     """The residual departures of the two measured ratios from four, $0.1\\%$ and $0.4\\%$, are not\ninterpreted.""", 1),
    ("""$1.2\\times10^{-10}$ and $A$ invariant to better than $10^{-10}$ relative --- a\nsolver tolerance, not an independent nonlinear measurement of $\\Rey_c$.""",
     """$1.2\\times10^{-10}$ and $A$ invariant to better than $10^{-10}$ relative, which is a\nsolver tolerance and not an independent nonlinear measurement of $\\Rey_c$.""", 1),

    # --- "exact" reserved for the theorem: measured agreements are given as numbers
    ("""ratio of $q_+-q_-$ stays within $3.7\\times10^{-4}$ of the exact value.""",
     """ratio of $q_+-q_-$ stays within $3.7\\times10^{-4}$ of the predicted value.""", 1),
    ("""Figure~\\ref{fig:sin_phi} shows the exact phase selection and""",
     """Figure~\\ref{fig:sin_phi} shows the phase selection predicted by \\eqref{eq:sin_phi} and""", 1),
    ("""all return the predicted reversal sign, the Grashof ratio exact to""",
     """all return the predicted reversal sign, the Grashof ratio agreeing to""", 1),
    ("""as consistent with exact Grashof linearity at the tolerance the runs reach, and""",
     """as consistent with Grashof linearity at the tolerance the runs reach, and""", 1),
    ("""artifact of the perturbative truncation. The exact reflection properties survive""",
     """artifact of the perturbative truncation. The reflection properties forced by symmetry survive""", 1),
]

edits_supp = [
    ("""The free-intercept fit is effectively exact and the response vanishes in the forced-convection limit""",
     """The free-intercept fit leaves an intercept of $-5.51\\times10^{-16}$ and the response vanishes in the forced-convection limit""", 1),
    ("""solver. Each curve crosses zero exactly once.""",
     """solver. Each curve crosses zero once.""", 1),
    ("""This distinction is important because the buoyancy term in the linearized
momentum equation is proportional to $\\Theta'(Y)$, rather than to the
temperature amplitude itself.""",
     """The buoyancy term in the linearized momentum equation is proportional to
$\\Theta'(Y)$ and not to the temperature amplitude itself.""", 1),
    ("""nonlinear mechanism redistributes the flux rather than creating an""",
     """nonlinear mechanism redistributes the flux; it creates no""", 1),
]


def main():
    for name, edits in ((MAIN, edits_main), (SUPP, edits_supp)):
        print("--- %s" % name)
        text = load(name)
        for old, new, count in edits:
            text = apply(text, old, new, name, count)
        io.open(os.path.join(ROOT, name), "w", encoding="utf-8",
                newline="\n").write(text)
        print("  written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

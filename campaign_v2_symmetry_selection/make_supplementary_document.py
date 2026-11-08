#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Split main_v2.tex into a main manuscript and a supplementary document.

Mechanical part only: it moves whole figure/table environments and the appendix
block out of main_v2.tex into supplementary_text.tex, and prints where each block
sat so the pointer sentences can be written by hand afterwards. Re-running is
idempotent from the current main_v2.tex: anything already absent is skipped.

    py -X utf8 campaign_v2_symmetry_selection/make_supplementary_document.py
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MAIN = os.path.join(ROOT, "main_v2.tex")
SUPP = os.path.join(ROOT, "supplementary_text.tex")

# Floats judged secondary for the main narrative; each stays cited from the text.
MOVE_FIGURES = [
    "fig:logic_diagram",
    "fig:Gr_sweep",
    "fig:dns_amplitude",
    "fig:dns_spatial",
    "fig:dns_timestep",
    "fig:dns_gr",
    "fig:dns_re",
    "fig:Rec_kartsch",
    "fig:Gr_invariance",
    "fig:phase_crossing",
    "fig:physical_temperature",
    "fig:physical_hydrodynamics",
    "fig:symmetry_forcing",
    "fig:wall_flux_physical",
]
MOVE_TABLES = [
    "tab:fixed_pe_Re",
    "tab:small_k_fit",
    "tab:dns_re_sweep",
    "tab:Gr_invariance",
    "tab:term_balance",
    "tab:paper_comparison",
]

# Groups in the supplement: (section title, [labels in presentation order])
GROUPS = [
    ("Perturbative parameter maps and scaling checks",
     ["tab:fixed_pe_Re", "tab:small_k_fit", "fig:Gr_sweep",
      "fig:Rec_kartsch", "fig:Gr_invariance", "tab:Gr_invariance"]),
    ("Verification of the nonlinear DNS",
     ["fig:dns_amplitude", "fig:dns_spatial", "fig:dns_timestep",
      "fig:dns_gr", "fig:dns_re", "tab:dns_re_sweep"]),
    ("Phase decomposition, term balance and transport comparison",
     ["fig:logic_diagram", "fig:phase_crossing", "tab:term_balance",
      "tab:paper_comparison"]),
    ("Two-dimensional physical structure of the mechanism",
     ["fig:physical_temperature", "fig:physical_hydrodynamics",
      "fig:symmetry_forcing", "fig:wall_flux_physical"]),
]

PREAMBLE = r"""\documentclass[12pt,a4paper]{article}

\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage[margin=2.5cm]{geometry}
\usepackage{setspace}
\usepackage{natbib}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{float}
\usepackage{siunitx}
\usepackage{tikz}
\usetikzlibrary{arrows.meta, positioning}
\usepackage[expansion=false]{microtype}   % expansion disabled: MiKTeX bitmap font fallback
\usepackage{xr}

\captionsetup{justification=raggedright, singlelinecheck=false}

\onehalfspacing
\bibliographystyle{apalike}   % single declaration

\newcommand{\Rey}{\mathit{Re}}
\newcommand{\Pe}{\mathit{Pe}}
\newcommand{\Grr}{\mathit{Gr}}
\newcommand{\Nu}{\mathit{Nu}}
\renewcommand{\Pr}{\mathit{Pr}}
\newcommand{\PiT}{\Pi_T}

\newtheorem{theorem}{Theorem}
\newtheorem{proposition}{Proposition}
\newtheorem{corollary}{Corollary}
\newtheorem{lemma}{Lemma}
\newenvironment{acknowledgments}{\section*{Acknowledgments}}{}
\newenvironment{ruledtabular}{}{}

% Supplementary floats and sections carry an explicit S so that a citation from
% the main manuscript can never be confused with a main-text number. The
% cross-references resolve in both directions once both documents have been
% compiled; the guard keeps a single-document build from failing.
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\IfFileExists{build_v2/main_v2.aux}
  {\externaldocument{build_v2/main_v2}}{}

\begin{document}

\title{\textbf{Supplementary Material for ``Symmetry-Selected Nonlinear Thermal
       Pumping in Periodically Heated Mixed-Convection Channels''}}
\author{BENIAICHE Ahmed \\
\textit{Laboratory of Fluid Mechanics, Ecole Militaire Polytechnique, Algiers, Algeria} \\
\textit{Email: beniaiche\_ahmed25@yahoo.fr}}
\date{}
\maketitle

\noindent This document collects the derivations, the verification records and the
supporting figures and tables of the study. It is referenced from the main
manuscript throughout: every item below is cited there, and its numbering carries
an \textbf{S} (\mbox{Fig.~S1}, \mbox{Table~S1}, Sec.~S1) so that the two documents
stay distinguishable in the text. Equation, section and theorem numbers quoted
inside are those of the main manuscript.
"""


def float_block(src, label, kind):
    """Return (start, end) character offsets of the environment holding label."""
    env = "figure" if kind == "fig" else "table"
    for m in re.finditer(r"\\begin\{%s\}" % env, src):
        depth = 1
        pos = m.end()
        while depth:
            nxt_b = src.find("\\begin{%s}" % env, pos)
            nxt_e = src.find("\\end{%s}" % env, pos)
            if nxt_e < 0:
                raise RuntimeError("unbalanced %s for %s" % (env, label))
            if 0 <= nxt_b < nxt_e:
                depth += 1
                pos = nxt_b + len(env) + 7
            else:
                depth -= 1
                pos = nxt_e + len(env) + 6
        start = src.rfind("\n", 0, m.start()) + 1
        end = src.find("\n", pos - 1)
        end = len(src) if end < 0 else end + 1
        if "\\label{%s}" % label in src[m.start():pos]:
            return start, end
    raise RuntimeError("no %s block for %s" % (env, label))


def main():
    src = io.open(MAIN, encoding="utf-8").read()
    blocks = {}

    # --- appendix block: from \appendix to the bibliography comment -------
    ap = src.find("\n\\appendix\n")
    if ap < 0:
        sys.exit("no \\appendix in the current main: appendices already moved")
    bib = src.find("%--------------------------------------------------------------\n% Bibliography",
                   ap)
    if bib < 0:
        sys.exit("cannot locate the bibliography marker")
    appendix = src[ap:bib]
    blocks["appendix"] = appendix
    print("appendix block: %d lines" % appendix.count("\n"))

    # --- floats -------------------------------------------------------------
    for kind, names in (("fig", MOVE_FIGURES), ("tab", MOVE_TABLES)):
        for label in names:
            try:
                s, e = float_block(src, label, kind)
            except RuntimeError as err:
                print("skip  %-28s (%s)" % (label, err))
                continue
            line = src.count("\n", 0, s) + 1
            blocks[label] = src[s:e]
            print("move  %-28s at main line %d (%d lines)"
                  % (label, line, src[s:e].count("\n")))

    # --- rewrite the main ---------------------------------------------------
    out = src
    # strip the appendix together with its surrounding comment rules
    cut = src[ap:bib]
    marker = "%--------------------------------------------------------------\n"
    tail_start = out.find(cut)
    out = out[:tail_start] + out[tail_start + len(cut):]
    for key, text in blocks.items():
        if key == "appendix":
            continue
        i = out.find(text)
        if i < 0:
            sys.exit("block %s disappeared while cutting" % key)
        out = out[:i] + out[i + len(text):]
    # collapse the blank lines left behind
    out = re.sub(r"\n{3,}", "\n\n", out)
    io.open(MAIN, "w", encoding="utf-8", newline="\n").write(out)
    print("main rewritten: %d lines (was %d)"
          % (out.count("\n"), src.count("\n")))

    # --- build the supplement ----------------------------------------------
    parts = [PREAMBLE]
    for title, labels in GROUPS:
        parts.append("\n\\section{%s}\n" % title)
        for label in labels:
            if label not in blocks:
                continue
            parts.append(blocks[label].rstrip() + "\n")
    app = blocks["appendix"].replace("\n\\appendix\n", "\n")
    app = app.replace(marker, "")
    parts.append("\n\\section{Derivations and independent checks}\n")
    parts.append(re.sub(r"\\section\{", r"\\subsection{", app).strip() + "\n")
    parts.append("""
%--------------------------------------------------------------
% Bibliography (single style already declared in the preamble)
\\bibliography{references}

\\end{document}
""")
    io.open(SUPP, "w", encoding="utf-8", newline="\n").write("".join(parts))
    print("supplementary_text.tex written: %d lines"
          % "".join(parts).count("\n"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Second split pass: move the two long exposition blocks out of main_v2.tex.

Cuts a block delimited by a start marker line and an end marker line (exclusive)
and appends it to the supplement under a section title. Run before writing the
pointer sentences, then check the diff.

    py -X utf8 campaign_v2_symmetry_selection/move_prose_sections.py
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MAIN = os.path.join(ROOT, "main_v2.tex")
SUPP = os.path.join(ROOT, "supplementary_text.tex")

BLOCKS = [
    ("\\subsection{Two independent mixed kernels}",
     "\\subsection{Numerical values and internal validation}",
     "Construction of the two independent mixed kernels",
     "sec:kernels"),
    ("\\subsection{Physical realization of the symmetry-selected pumping mechanism}",
     "\\subsection{Buoyancy as the hydrodynamic mediator}",
     "Physical realization of the symmetry-selected mechanism",
     "sec:physical_mechanism"),
]


def cut(src, start_marker, end_marker):
    i = src.find(start_marker)
    if i < 0:
        return None, src, None
    j = src.find(end_marker, i)
    if j < 0:
        sys.exit("end marker %s not found after %s" % (end_marker, start_marker))
    block = src[i:j]
    # keep the anchor of the following heading stable
    return block, src[:i] + src[j:], i


def main():
    src = io.open(MAIN, encoding="utf-8").read()
    supp = io.open(SUPP, encoding="utf-8").read()
    moved = []
    for start_marker, end_marker, title, anchor in BLOCKS:
        block, src, pos = cut(src, start_marker, end_marker)
        if block is None:
            print("skip  %s (already moved)" % title)
            continue
        block = re.sub(r"\\subsection\{", "\\\\section{", block, count=1)
        block = block.replace("\\subsubsection{", "\\subsection{")
        if "\\label{%s}" % anchor not in block:
            block = "\\label{%s}\n" % anchor + block
        block = block.rstrip() + "\n"
        moved.append((title, block, pos))
        print("move  %-58s at main line %d (%d lines)"
              % (title, src[:pos].count("\n") + 1, block.count("\n")))

    if not moved:
        print("nothing to do")
        return 0

    insert_at = supp.find("\n\\bibliography{references}")
    if insert_at < 0:
        insert_at = supp.find("%--------------------------------------------------------------\n% Bibliography")
    if insert_at < 0:
        sys.exit("cannot locate the bibliography block of the supplement")
    addition = ""
    for title, block, _ in moved:
        addition += "\n%% ---- %s ----\n" % title + block
    supp = supp[:insert_at] + addition + supp[insert_at:]

    io.open(MAIN, "w", encoding="utf-8", newline="\n").write(
        re.sub(r"\n{3,}", "\n\n", src))
    io.open(SUPP, "w", encoding="utf-8", newline="\n").write(supp)
    print("main now %d lines, supplement %d lines"
          % (src.count("\n"), supp.count("\n")))
    return 0


if __name__ == "__main__":
    sys.exit(main())

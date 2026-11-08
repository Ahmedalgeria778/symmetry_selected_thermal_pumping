"""Notation consistency audit for main_v2.tex (plan item B7).

Three documented conventions are enforced across the manuscript:
  * dimensionless groups are the macros \\Rey, \\Pe, \\Grr, \\Pr -- never raw
    italic Re/Pe/Gr/Pr inside maths;
  * the Fourier reconstruction is Re[. e^{-ikX}] everywhere;
  * Pi_T names the linear pumping coefficient and Nu_{2,0} the wall-flux
    response; the two are never equated without the amplitude factors.
Exit code is non-zero when a convention is violated.
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.abspath(os.path.join(HERE, os.pardir, "main_v2.tex"))

GROUPS = {"Re": r"\Rey", "Pe": r"\Pe", "Gr": r"\Grr", "Pr": r"\Pr"}
RAW = re.compile(r"(?<![A-Za-z_\\])(Pe|Re|Gr|Pr)(?![A-Za-z_])")


def strip_comments(text):
    return "\n".join(re.sub(r"(?<!\\)%.*", "", line) for line in text.split("\n"))


def math_spans(text):
    spans = []
    for m in re.finditer(r"\$[^$]*\$", text):
        spans.append((m.start(), m.group()))
    for env in ("equation", "align", "gather", "array", "tabular", "tabular*"):
        for m in re.finditer(r"\\begin\{%s\*?\}.*?\\end\{%s\*?\}" % (env, env),
                             text, re.S):
            spans.append((m.start(), m.group()))
    return spans


def main():
    text = strip_comments(io.open(TEX, encoding="utf-8").read())

    def lineno(pos):
        return text.count("\n", 0, pos) + 1

    failures = []

    # 1. raw italic dimensionless groups inside maths
    for start, chunk in math_spans(text):
        chunk = re.sub(r"\\operatorname\{Re\}|\\Re\b|\\text\{[^}]*\}", " ", chunk)
        for m in RAW.finditer(chunk):
            failures.append((lineno(start), "raw italic %s; use %s"
                             % (m.group(0), GROUPS[m.group(0)][1:])))

    # 2. Fourier reconstruction sign convention
    for m in re.finditer(r"e\^\{ikX\}", text):
        failures.append((lineno(m.start()), "e^{ikX} reconstruction; the paper "
                         "fixes e^{-ikX}"))

    # 3. coefficient vs response: Nu_{2,0} = Pi_T with no amplitude factors
    for m in re.finditer(r"\\Nu_\{?2,0\}?\s*(\(\S{1,12}\))?\s*=\s*\\PiT(?!\^)", text):
        failures.append((lineno(m.start()),
                             "Nu_2,0 equated to Pi_T without amplitude factors"))

    # 4. tables: report which pumping symbol each table speaks about
    print("table inventory")
    for m in re.finditer(r"\\begin\{table\}.*?\\end\{table\}", text, re.S):
        blk = m.group(0)
        lab = re.search(r"\\label\{(tab:[^}]*)\}", blk)
        pi = len(re.findall(r"\\PiT|\\Pi_\{?T", blk))
        nu = len(re.findall(r"\\Nu_\{?2,0", blk))
        a = len(re.findall(r"\bA=_?", blk))
        print("  %-28s Pi_T:%-3d Nu_2,0:%-3d A:%d"
              % (lab.group(1) if lab else "(unlabelled)", pi, nu, a))

    print("\nreal-part operator usage")
    print("  \\operatorname{Re}: %d" % len(re.findall(r"\\operatorname\{Re\}", text)))
    print("  \\Re              : %d" % len(re.findall(r"\\Re(?![A-Za-z])", text)))

    print("\ndeclared macros actually used")
    for name in ("Pec", "Gra"):
        n = len(re.findall(r"\\" + name + r"(?![A-Za-z])", text))
        # one hit is the \newcommand definition itself
        print("  \\%-4s occurrences: %d %s" % (name, n,
              "(definition only -- dead)" if n <= 1 else ""))

    if failures:
        print("\nFAIL: %d notation violations" % len(failures))
        for ln, msg in failures:
            print("  line %4d  %s" % (ln, msg))
        return 1
    print("\nPASS: notation conventions consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())

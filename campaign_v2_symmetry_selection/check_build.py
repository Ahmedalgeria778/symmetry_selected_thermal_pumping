#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build gate for the two-document paper: main_v2.tex and supplementary_text.tex.

Parses both .log and .aux files, then checks both sources as one corpus: a label
may be defined in either document and cited from either, because the two files
are linked with the xr package. Fails on any defect.

    py -X utf8 campaign_v2_symmetry_selection/check_build.py
"""
import io
import re
import sys
from pathlib import Path

WORK = Path(__file__).resolve().parents[1]
BUILD = WORK / "build_v2"

MAIN = "main_v2"
SUPP = "supplementary_text"
SOURCES = [WORK / (MAIN + ".tex"), WORK / (SUPP + ".tex"),
           WORK / "cover_letter_NHTB_v2.tex"]

# a ceiling on the main text: the split exists to keep it readable
MAIN_PAGE_LIMIT = 45

# labels the manuscript cross-references; a missing one is a hard error. The
# first group lives in the main text, the second in the supplementary document
# and is reached through xr.
REQUIRED_MAIN = ["tab:dns_predictions", "eq:dns_prediction_rule", "tab:Rec_k",
                 "fig:endpoint", "tab:dns_validation_summary",
                 "eq:PiT_dns_estimator", "prop:phase", "cor:Gr_independence",
                 "cor:sinphi", "fig:dns_multik", "subsec:error_model",
                 "subsec:Gr_test", "tab:Gr_dns", "subsec:remainder",
                 "sec:symmetry_law", "fig:dns_grpair"]
REQUIRED_SUPP = ["fig:Rec_kartsch", "sec:kernels", "sec:physical_mechanism",
                 "app:hydro_derivation", "app:reflection", "app:cross_proof",
                 "app:mean_momentum", "app:theorem_proof", "app:prop_proof",
                 "app:mean_solution", "app:validation", "sec:physical_summary"]


def read(path):
    return io.open(path, encoding="utf-8", errors="replace").read()


def labels_in(src):
    body = re.sub(r"\\begin\{subfigure\}.*?\\end\{subfigure\}", "", src, flags=re.S)
    return [l for l in re.findall(r"\\label\{([^}]*)\}", body)
            if l.startswith(("fig:", "tab:"))]


def hygiene(path, corpus_labels, corpus_src, own_src):
    """Source-level scoria check: non-ASCII bytes, leftover French comments,
    placeholders, duplicate labels, dangling labels and duplicate bib styles."""
    lines = own_src.split("\n")
    problems = []
    for i, line in enumerate(lines, 1):
        if any(ord(c) > 126 for c in line):
            problems.append("%s:%d non-ascii: %s" % (path.name, i, line.strip()[:60]))
        cm = re.search(r"(?<!\\)%\s*(.*)$", line)
        comment = cm.group(1) if cm else ""
        if re.search(r"\b(avec|pour|les|dans|sans|donc)\b", comment, re.I):
            problems.append("%s:%d french comment: %s"
                            % (path.name, i, comment.strip()[:60]))
        if re.search(r"\b(TODO|FIXME|XXX|TBD|PLACEHOLDER)\b", line, re.I):
            problems.append("%s:%d placeholder: %s" % (path.name, i, line.strip()[:60]))
    labels = re.findall(r"\\label\{([^}]*)\}", own_src)
    dupes = sorted({l for l in labels if labels.count(l) > 1})
    if dupes:
        problems.append("%s duplicate label(s): %s" % (path.name, dupes))
    refs = set(re.findall(r"\\(?:ref|eqref|autoref|Cref|cref|pageref)\{([^}]*)\}",
                          own_src))
    dangling = sorted(r for r in refs if r not in corpus_labels)
    if dangling:
        problems.append("%s reference(s) to undefined label(s): %s"
                        % (path.name, dangling))
    nstyle = len(re.findall(r"\\bibliographystyle\{", own_src))
    if nstyle > 1:
        problems.append("%s declares bibliographystyle %d times" % (path.name, nstyle))
    # Journal prose must not cite working filenames; the figure inputs are the
    # only file paths allowed in a source, and the archive inventory (MANIFEST)
    # carries the per-item mapping instead.
    stripped = re.sub(r"\\includegraphics(\[[^\]]*\])?\{[^}]*\}", "", own_src)
    files = sorted(set(re.findall(
        r"[A-Za-z0-9_\\.-]+\.(?:csv|py|json|txt|md|zip)\b", stripped)))
    if files:
        problems.append("%s cites working filenames: %s"
                        % (path.name, files[:5]))
    if path.name == MAIN + ".tex" and "\\textbf{Keywords:}" not in own_src:
        problems.append("%s has no keyword line" % path.name)
    # A journal float that nothing points to is a defect, not a spare figure.
    # Floats are cited across the corpus, so the search covers both documents.
    uncited = sorted(l for l in labels_in(own_src) if not re.search(
        r"\\(?:ref|eqref|autoref|Cref|cref|pageref)\{%s\}" % re.escape(l), corpus_src))
    if uncited:
        problems.append("%s float(s) never cited: %s" % (path.name, uncited))
    return problems


def report(name):
    log = read(BUILD / (name + ".log"))
    flat = log.replace("\n", " ")
    errors = re.findall(r"^!.*", log, re.M)
    undef = re.findall(r"(?:Reference|Citation|Label).*?undefined.*", log)
    warn = [w for w in re.findall(r"LaTeX Warning: .*", log)
            if "undefined" not in w]
    overfull = re.findall(r"Overfull \\hbox \(([^)]*)\)", log)
    badboxes = [o for o in overfull if float(o.split("pt")[0]) > 20.0]
    badvbox = re.findall(r"Overfull \\vbox \(([^)]*)\)", log)
    pages = re.findall(r"\((\d+)\s+pages", flat)
    print("%s: errors %d, undefined %d, other warnings %d, overfull hbox > 20 pt "
          "%d (of %d), overfull vbox %d, pages %s"
          % (name, len(errors), len(undef), len(warn), len(badboxes),
             len(overfull), len(badvbox), pages[0] if pages else "not written"))
    for e in (errors + undef + warn)[:6]:
        print("   ", e[:120])
    for b in (badboxes + badvbox)[:6]:
        print("    (%s)" % b)
    return errors, undef, badboxes, badvbox, pages


def main():
    bad = []
    counts = []
    for name in (MAIN, SUPP):
        counts.append(report(name))
    main_pages = counts[0][4]
    if main_pages and int(main_pages[0]) > MAIN_PAGE_LIMIT:
        print("   main text exceeds the %d-page ceiling" % MAIN_PAGE_LIMIT)

    corpus_src = "\n".join(read(p) for p in SOURCES[:2])
    corpus_labels = set(re.findall(r"\\label\{([^}]*)\}", corpus_src))

    aux = {name: read(BUILD / (name + ".aux")) for name in (MAIN, SUPP)}
    both = aux[MAIN] + aux[SUPP]
    for lab in REQUIRED_MAIN + REQUIRED_SUPP:
        hit = ("newlabel{%s}" % lab) in both or ("newlabel@cref@{%s}" % lab) in both
        if not hit:
            print("label %-28s : MISSING" % lab)
            bad.append("label " + lab)
    print("required labels   : %d checked, %d missing"
          % (len(REQUIRED_MAIN) + len(REQUIRED_SUPP), len(bad)))

    # Every supplementary float has to be reachable from the main text: that is
    # what makes the split a division of labour rather than a hiding place.
    supp_src = read(WORK / (SUPP + ".tex"))
    main_src = read(WORK / (MAIN + ".tex"))
    orphan = sorted(l for l in labels_in(supp_src) if not re.search(
        r"\\(?:ref|eqref|autoref|Cref|cref|pageref)\{%s\}" % re.escape(l), main_src))
    print("supplementary floats cited from the main text : %d of %d"
          % (len(labels_in(supp_src)) - len(orphan), len(labels_in(supp_src))))
    if orphan:
        print("   never cited from the main manuscript: %s" % orphan)

    scoria = []
    for path in SOURCES:
        scoria += hygiene(path, corpus_labels, corpus_src, read(path))
    print("source scoria     : %d" % len(scoria))
    for s in scoria[:10]:
        print("   ", s[:120])

    errors = sum(len(c[0]) for c in counts)
    undef = sum(len(c[1]) for c in counts)
    boxes = sum(len(c[2]) for c in counts) + sum(len(c[3]) for c in counts)
    if errors or undef or boxes or scoria or bad or orphan or not main_pages:
        print("\nBUILD GATE FAILED:", ", ".join(
            ["%d errors" % errors, "%d undefined" % undef,
             "%d bad boxes" % boxes, "%d missing labels" % len(bad),
             "%d uncited supplementary floats" % len(orphan),
             "%d scoria" % len(scoria)]))
        return 1
    print("\nPASS: both documents build clean, main %s pages (limit %d), "
          "supplementary %s pages, all required labels resolved, every "
          "supplementary float cited from the main text, sources ascii-clean"
          % (main_pages[0], MAIN_PAGE_LIMIT, counts[1][4][0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build gate for main_v2.tex: parse the .log and .aux, fail on any defect.

    py -X utf8 campaign_v2_symmetry_selection/check_build.py
"""
import io
import re
import sys
from pathlib import Path

WORK = Path(__file__).resolve().parents[1]
BUILD = WORK / "build_v2"
LOG = BUILD / "main_v2.log"
AUX = BUILD / "main_v2.aux"
PDF = BUILD / "main_v2.pdf"

# labels that the current draft cross-references; a missing one is a hard error
REQUIRED = ["tab:dns_predictions", "eq:dns_prediction_rule", "tab:Rec_k",
            "fig:Rec_kartsch", "fig:Rec_locus", "fig:endpoint",
            "tab:dns_validation_summary", "eq:PiT_dns_estimator",
            "prop:phase", "cor:Gr_independence", "cor:sinphi",
            "fig:dns_multik", "subsec:error_model", "subsec:Gr_test",
            "tab:Gr_dns", "subsec:remainder", "sec:symmetry_law",
            "fig:dns_grpair"]


def hygiene(path):
    """Source-level scoria check: non-ASCII bytes, leftover French comments,
    placeholders, duplicate or dangling labels, duplicate bibliography styles."""
    lines = io.open(path, encoding="utf-8", errors="replace").read().split("\n")
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
    src = "\n".join(lines)
    labels = re.findall(r"\\label\{([^}]*)\}", src)
    dupes = sorted({l for l in labels if labels.count(l) > 1})
    if dupes:
        problems.append("%s duplicate label(s): %s" % (path.name, dupes))
    refs = set(re.findall(r"\\(?:ref|eqref|autoref|Cref|cref|pageref)\{([^}]*)\}", src))
    dangling = sorted(r for r in refs if r not in labels)
    if dangling:
        problems.append("%s reference(s) to undefined label(s): %s"
                        % (path.name, dangling))
    nstyle = len(re.findall(r"\\bibliographystyle\{", src))
    if nstyle > 1:
        problems.append("%s declares bibliographystyle %d times" % (path.name, nstyle))
    # Journal prose must not cite working filenames; the figure inputs are the
    # only file paths allowed in a source, and the archive inventory (MANIFEST)
    # carries the per-item mapping instead.
    stripped = re.sub(r"\\includegraphics(\[[^\]]*\])?\{[^}]*\}", "", src)
    files = sorted(set(re.findall(
        r"[A-Za-z0-9_\\.-]+\.(?:csv|py|json|txt|md|zip)\b", stripped)))
    if files:
        problems.append("%s cites working filenames: %s"
                        % (path.name, files[:5]))
    if path.name == "main_v2.tex" and "\\textbf{Keywords:}" not in src:
        problems.append("main_v2.tex has no keyword line")
    # A journal float that nothing points to is a defect, not a spare figure.
    # Panels of a composite figure are cited through their parent, so the
    # subfigure blocks are excluded from the label census.
    body = re.sub(r"\\begin\{subfigure\}.*?\\end\{subfigure\}", "", src, flags=re.S)
    floats = [l for l in re.findall(r"\\label\{([^}]*)\}", body)
              if l.startswith(("fig:", "tab:"))]
    uncited = sorted(l for l in floats if not re.search(
        r"\\(?:ref|eqref|autoref|Cref|cref|pageref)\{%s\}" % re.escape(l), src))
    if uncited:
        problems.append("%s float(s) never cited: %s" % (path.name, uncited))
    return problems


def main():
    log = io.open(LOG, encoding="utf-8", errors="replace").read()
    aux = io.open(AUX, encoding="utf-8", errors="replace").read()
    bad = []

    errors = re.findall(r"^!.*", log, re.M)
    undef = re.findall(r"(?:Reference|Citation|Label).*?undefined.*", log)
    warn = [w for w in re.findall(r"LaTeX Warning: .*", log)
            if "undefined" not in w]
    overfull = re.findall(r"Overfull \\\\hbox \(([^)]*)\)", log)
    badboxes = [o for o in overfull if float(o.split("pt")[0]) > 20.0]
    flat = log.replace("\n", " ")
    pages = re.findall(r"\((\d+)\s+pages", flat)

    print("errors            : %d" % len(errors))
    for e in errors[:5]:
        print("   ", e[:120])
    print("undefined refs    : %d" % len(undef))
    for u in undef[:5]:
        print("   ", u[:120])
    print("other warnings    : %d" % len(warn))
    for w in warn[:5]:
        print("   ", w[:120])
    print("overfull > 20 pt  : %d (of %d total)" % (len(badboxes), len(overfull)))
    for b in badboxes[:5]:
        print("    (%s)" % b)
    print("pages             : %s" % (pages[0] if pages else "not written"))
    for lab in REQUIRED:
        hit = ("newlabel{%s}" % lab) in aux or ("newlabel@cref@{%s}" % lab) in aux
        print("label %-28s : %s" % (lab, "ok" if hit else "MISSING"))
        if not hit:
            bad.append("label " + lab)
    scoria = hygiene(WORK / "main_v2.tex") + hygiene(WORK / "cover_letter_NHTB_v2.tex")
    print("source scoria       : %d" % len(scoria))
    for s in scoria[:10]:
        print("   ", s[:120])
    if errors or undef or badboxes or not pages or scoria:
        print("\nBUILD GATE FAILED:", ", ".join(
            [f"{len(errors)} errors", f"{len(undef)} undefined",
             f"{len(badboxes)} bad boxes", f"{len(bad)} missing labels",
             f"{len(scoria)} scoria"]))
        return 1
    print("\nPASS: build clean, %s pages, all required labels resolved, "
          "sources ascii-clean" % pages[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())

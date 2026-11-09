"""One-off check of the compiled cover letter: log parse and page extraction.

The letter is the one-page document requested by the journal, so the gate
checks two things at once: the statements a submission letter has to carry
(NEEDLES) and the promotional framing the referee asked to remove (FORBIDDEN).
Both lists are quoted from the compiled text, so a silent edit fails here.
"""
import io
import re
import sys
from pathlib import Path

from pypdf import PdfReader

WORK = Path(__file__).resolve().parents[1]

MAX_PAGES = 1
NEEDLES = [
    # venue, article type and the title exactly as the manuscript spells it
    "Numerical Heat Transfer, Part B",
    "Taylor & Francis",
    "Research Article",
    "Mixed-Convection Channel",
    # the three physical statements the letter commits to
    "reverses with the Reynolds number",
    "adds a quadrature coefficient",
    "fully nonlinear stationary DNS",
    # the campaign result, stated without the promotional frame
    "twelve non-degenerate cases",
    "exceeds the error bound",
    # the companion disclosure and the two-document split
    "submitted to this journal on 1 November 2026",
    "no result, figure or derivation of that study is reproduced",
    "separately numbered supplementary document",
    "S-numbering",
    # the public author deposit, offered without asking the editor to rule on it
    "symmetry_selected_thermal_pumping",
    "figure-generation scripts",
    "available for the editor",
    # the routine declarations
    "not under consideration elsewhere",
    "no conflict of interest",
    "no external funding",
]
FORBIDDEN = [
    "headline", "pillars", "falsif", "destructive", "sharpest",
    "twelve-word", "64 files", "interleaved", "author-and-affiliation",
    "refined-streamwise-mesh",
]

log = io.open(WORK / "build_cover" / "cover_letter_NHTB_v2.log",
              encoding="utf-8", errors="replace").read()
flat = log.replace("\n", " ")
errors = re.findall(r"^!.*", log, re.M)
undef = re.findall(r"(?:Reference|Citation).*?undefined.*", log)
boxes = [float(o.split("pt")[0]) for o in
         re.findall(r"Overfull \\hbox \(([^)]*)\)", log)]
pages = re.findall(r"\((\d+)\s+pages?\b", flat)
print("errors=%d undefined=%d overfull>20pt=%d pages=%s"
      % (len(errors), len(undef), len([b for b in boxes if b > 20.0]),
         pages[0] if pages else "?"))
for e in errors[:5]:
    print("   ", e[:120])

pdf = PdfReader(str(WORK / "build_cover" / "cover_letter_NHTB_v2.pdf"))
text = re.sub(r"\s+", " ", " ".join((p.extract_text() or "") for p in pdf.pages))
# pdflatex breaks words across lines, so the extracted text contains
# "ex- ceeds" and "docu- ment"; needles are matched against the
# de-hyphenated variant as well, since a real hyphen is the only cost there.
joined = re.sub(r"([a-z])- ([a-z])", r"\1\2", text)


def found(needle):
    return needle in text or needle.lower() in joined.lower()


print("characters: %d" % len(text))
missing = [n for n in NEEDLES + list(sys.argv[1:]) if not found(n)]
present = [f for f in FORBIDDEN if f in text or f in joined]
for needle in NEEDLES:
    print("%-52s %s" % (needle[:52], "FOUND" if found(needle) else "MISSING"))
for phrase in present:
    print("%-52s FORBIDDEN" % phrase)

too_long = bool(pages) and int(pages[0]) > MAX_PAGES
problems = (list(errors) + list(undef)
            + ["overfull box %s pt" % b for b in boxes if b > 20.0]
            + ["missing statement: %s" % m for m in missing]
            + ["promotional wording left in the letter: %s" % p for p in present]
            + (["letter is %s pages, limit is %d" % (pages[0], MAX_PAGES)]
               if too_long else []))
if problems:
    print("\nCOVER LETTER GATE FAILED: %d defect(s)" % len(problems))
    for p in problems[:8]:
        print("   ", str(p)[:120])
    sys.exit(1)
print("\nPASS: cover letter compiles clean, %s page (limit %d), all quoted "
      "statements present, no promotional wording"
      % (pages[0] if pages else "?", MAX_PAGES))

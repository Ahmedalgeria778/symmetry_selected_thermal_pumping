"""One-off check of the compiled cover letter: log parse and page extraction."""
import io
import re
import sys
from pathlib import Path

from pypdf import PdfReader

WORK = Path(__file__).resolve().parents[1]

# Statements the letter must contain once the Gr-pair campaign and the
# falsification of the error model are integrated; they are quoted from the
# compiled text so a silent edit to the letter fails this gate. The last four
# entries are the venue framing and the study-1 disclosure for the NHT-B
# submission. The page limit is the one requested for the submission.
MAX_PAGES = 2
NEEDLES = [
    "three mutually independent pillars",
    "headline number",
    "eleven digits",
    "4.0003",
    "4.0015",
    "3.996",
    "4.015",
    "Richardson",
    "falsified",
    "author-and-affiliation",
    "refined-streamwise-mesh",
    "three shared twelve-word shingles",
    "64 files",
    "Numerical Heat Transfer",
    "Taylor & Francis",
    "submitted to this journal on 1 November 2026",
    "shooting",
    # the public author deposit and the preprint-policy question it raises
    "symmetry_selected_thermal_pumping",
    "preprint policy",
]
log = io.open(WORK / "build_cover" / "cover_letter_NHTB_v2.log",
              encoding="utf-8", errors="replace").read()
flat = log.replace("\n", " ")
errors = re.findall(r"^!.*", log, re.M)
undef = re.findall(r"(?:Reference|Citation).*?undefined.*", log)
boxes = [float(o.split("pt")[0]) for o in
         re.findall(r"Overfull \\hbox \(([^)]*)\)", log)]
pages = re.findall(r"\((\d+)\s+pages", flat)
print("errors=%d undefined=%d overfull>20pt=%d pages=%s"
      % (len(errors), len(undef), len([b for b in boxes if b > 20.0]),
         pages[0] if pages else "?"))
for e in errors[:5]:
    print("   ", e[:120])

pdf = PdfReader(str(WORK / "build_cover" / "cover_letter_NHTB_v2.pdf"))
text = re.sub(r"\s+", " ", " ".join((p.extract_text() or "") for p in pdf.pages))
print("characters: %d" % len(text))
needles = NEEDLES + list(sys.argv[1:])
missing = [n for n in needles if n not in text]
for needle in needles:
    print("%-52s %s" % (needle[:52], "FOUND" if needle in text else "MISSING"))

too_long = bool(pages) and int(pages[0]) > MAX_PAGES
problems = (list(errors) + list(undef)
            + ["overfull box %s pt" % b for b in boxes if b > 20.0]
            + ["missing statement: %s" % m for m in missing]
            + (["letter is %s pages, limit is %d" % (pages[0], MAX_PAGES)]
               if too_long else []))
if problems:
    print("\nCOVER LETTER GATE FAILED: %d defect(s)" % len(problems))
    for p in problems[:8]:
        print("   ", str(p)[:120])
    sys.exit(1)
print("\nPASS: cover letter compiles clean, %s pages (limit %d), all quoted "
      "statements present" % (pages[0] if pages else "?", MAX_PAGES))

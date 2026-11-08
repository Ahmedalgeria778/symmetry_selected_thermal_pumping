"""Gate: the abstract obeys the journal rules and the end matter is complete.

Checks, on the manuscript source: the abstract is one paragraph of at most 250
words with no displayed equation, no citation and no cross-reference; the
document ends with a Supplementary Material section whose title is exactly that
phrase, followed by a disclosure statement, an author-contribution statement and
a funding statement; and the body points at the supplementary document with the
words "supplementary material".
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "main_v2.tex"
LIMIT = 250

text = MAIN.read_text(encoding="utf-8")
fail = []

abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", text, re.S)
if abstract is None:
    fail.append("no abstract environment")
else:
    body = abstract.group(1)
    if re.search(r"(?<!\S)\n\s*\n(?!\S)", "\n" + body + "\n"):
        fail.append("abstract is not a single paragraph")
    for name, pat in (("display math \\[", r"\\\["),
                      ("equation environment", r"\\begin\{equation\}"),
                      ("align environment", r"\\begin\{align"),
                      ("citation", r"\\cite"),
                      ("cross-reference", r"\\ref")):
        if re.search(pat, body):
            fail.append("abstract contains a " + name)
    words = re.findall(r"[A-Za-z0-9][A-Za-z0-9%\-'_.]*",
                       re.sub(r"\\[a-zA-Z]+\*?", " ", body))
    print("abstract words: %d (limit %d)" % (len(words), LIMIT))
    if len(words) > LIMIT:
        fail.append("abstract has %d words" % len(words))

for title in ("Supplementary Material", "Disclosure statement",
              "Author contributions", "Funding"):
    if ("\\section*{%s}" % title) not in text:
        fail.append("missing \\section*{%s}" % title)

concl = text.find("\\section{Conclusions}")
supp = text.find("\\section*{Supplementary Material}")
if concl >= 0 and supp >= 0 and supp < concl:
    fail.append("Supplementary Material section precedes the Conclusions")

pointers = len(re.findall(r"of the supplementary material", text))
print("body pointers using the phrase: %d" % pointers)
if pointers < 5:
    fail.append("only %d pointers say 'supplementary material'" % pointers)

for line in text.splitlines():
    if any(ord(c) > 127 for c in line):
        fail.append("non-ascii source line: " + line[:60])
        break

if fail:
    print("FAIL:")
    for f in fail:
        print("  - " + f)
    sys.exit(1)
print("PASS: abstract and end matter obey the editorial rules")

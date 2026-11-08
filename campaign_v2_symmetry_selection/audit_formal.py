import re, json, sys
from pypdf import PdfReader

PDF = r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification\main.pdf"
TEX = r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification\main.tex"
OUT = r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification\campaign_v2_symmetry_selection\B_formal_audit.txt"

tex = open(TEX, encoding="utf-8", errors="replace").read()
r = PdfReader(PDF)
pages = [p.extract_text() or "" for p in r.pages]
full = "\n".join(pages)
lines = []
def P(s):
    lines.append(s)

P("PDF pages: %d" % len(pages))

# B3 duplicate subsection headings in tex
subs = re.findall(r"\\subsection\{([^}]*)\}", tex)
dup = [s for s in set(subs) if subs.count(s) > 1]
P("\n[B3] duplicate \\subsection titles in .tex: %s" % (dup or "NONE"))
P("[B3] occurrences 'Solvability and wall-flux identities' in .tex: %d" % tex.count("Solvability and wall-flux identities"))
P("[B3] occurrences in PDF text: %d" % full.count("Solvability and wall-"))

# B4 'Consider' doubled
P("\n[B4] 'Consider' at line start in .tex: %d" % len(re.findall(r"\nConsider\b", tex)))
P("[B4] 'Consider Consider' in PDF: %s" % ("Consider Consider" in full))

# B5 corrupted schematic formula
P("\n[B5] 'partial \\gamma' in .tex: %d" % tex.count(r"\partial \gamma"))
P("[B5] arrow chains in .tex: %d" % len(re.findall(r"\\to\s", tex)))
seg = re.findall(r"[^\n]*\\to[^\n]*", tex)
P("[B5] lines containing \\to :")
for s in seg[:12]:
    P("      " + s.strip()[:140])

# B1 '1 1 1 ... 152' artefact: look for suspicious repeated '1' runs in extracted text
runs = [i for i, t in enumerate(pages) if re.search(r"(^|\s)(1\s+){5,}", t)]
P("\n[B1] pages with long runs of '1': %s" % ([i + 1 for i in runs] or "NONE"))
if runs:
    for i in runs[:3]:
        m = re.search(r"((1\s+){5,}.{0,40})", pages[i])
        P("      page %d context: %r" % (i + 1, (m.group(0)[:80] if m else "")))

# B2 nomenclature empty units
P("\n[B2] '[]' empty-unit tokens in .tex: %d" % len(re.findall(r"\[\s*\]", tex)))
P("[B2] nomenclature rows: %d" % len(re.findall(r"^\$", tex, re.M)))

# B10 equation numbering: collect tags
nums = re.findall(r"\\label\{eq:([^}]*)\}", tex)
P("\n[B10] number of equation labels: %d" % len(nums))
pdfnums = sorted(set(int(x) for x in re.findall(r"\((\d{1,3})\)", full)))
P("[B10] PDF rendered equation numbers range: %s..%s, count=%d" % (pdfnums[:1], pdfnums[-1:], len(pdfnums)))
missing = [n for n in range(1, (pdfnums[-1] if pdfnums else 0) + 1) if n not in pdfnums]
P("[B10] gaps in rendered numbering (1..max): %s" % (missing or "NONE"))

# B8 Fig 8 caption: find captions
caps = re.findall(r"\\caption\{([^}]{0,200})", tex, re.S)
P("\n[B8] captions count: %d" % len(caps))
for i, c in enumerate(caps):
    P("   [%2d] %s" % (i + 1, " ".join(c.split())[:150]))

# B6 references rendering
refstart = max(i for i, t in enumerate(pages) if "REFERENCES" in t.upper() or "Bibliography" in t)
P("\n[B6] reference section starts on page %d" % (refstart + 1))
P("[B6] raw rendered reference text:")
P(pages[refstart][:900])

# B11 email
P("\n[B11] email occurrences in .tex: %s" % set(re.findall(r"[\w.\-]+@[\w.\-]+", tex)))

# B7 tables: Pi_T vs Nu_{2,0} notation inside tables
tbl = re.findall(r"\\begin\{table\}(.*?)\\end\{table\}", tex, re.S)
P("\n[B7] table environments: %d" % len(tbl))
for i, t in enumerate(tbl):
    head = " ".join(t.split())[:110]
    P("   table %d: %s | has PiT=%s has Nu20=%s" % (i + 1, head,
        r"\PiT" in t or "Pi_T" in t, "Nu_{2,0}" in t or r"Nu_{2,0}" in t))

open(OUT, "w", encoding="utf-8").write("\n".join(lines))
print("wrote", OUT, len(lines), "lines")

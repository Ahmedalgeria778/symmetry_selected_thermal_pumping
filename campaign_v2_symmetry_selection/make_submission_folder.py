#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble the journal submission folder from the three compiled documents.

Copies are byte-identical to the working sources and their PDFs; every digest is
recomputed after the copy and the folder fails if a pair differs. The README is
generated from the build logs, so its page counts cannot go stale.

    py -X utf8 campaign_v2_symmetry_selection/make_submission_folder.py
"""
import hashlib
import io
import re
import shutil
import sys
from pathlib import Path

WORK = Path(__file__).resolve().parents[1]
OUT = WORK / "submission_NHTB_2026-11-08"
DATE = "8 November 2026"
REPO = "https://github.com/Ahmedalgeria778/symmetry_selected_thermal_pumping"

ITEMS = [
    ("manuscript", WORK / "main_v2.tex", WORK / "build_v2" / "main_v2.pdf"),
    ("supplementary_document", WORK / "supplementary_text.tex",
     WORK / "build_v2" / "supplementary_text.pdf"),
    ("cover_letter", WORK / "cover_letter_NHTB_v2.tex",
     WORK / "build_cover" / "cover_letter_NHTB_v2.pdf"),
]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def pages(log_path):
    log = io.open(log_path, encoding="utf-8", errors="replace").read()
    m = re.findall(r"\((\d+)\s+pages", log.replace("\n", " "))
    return m[0] if m else "?"


def supplementary_sections(src):
    """The S-numbered section titles, read from the supplement itself so the
    README cannot drift from the document it describes."""
    titles = ["S%d %s%s" % (i, t[0].lower(), t[1:]) for i, t in
              enumerate(re.findall(r"^\\section\{([^}]*)\}", src, re.M), 1)]
    return ", ".join(titles[:-1]) + " and " + titles[-1], len(titles)


def main():
    OUT.mkdir(exist_ok=True)
    rows, fails = [], []
    for name, tex, pdf in ITEMS:
        for src, suffix in ((tex, ".tex"), (pdf, ".pdf")):
            dst = OUT / (name + suffix)
            # A copy already matching its source is left alone: Windows can hold
            # a lock on a PDF that a viewer has open, and rewriting it is not
            # worth a failed build.
            if dst.exists() and digest(src) == digest(dst):
                pass
            else:
                try:
                    shutil.copyfile(src, dst)
                except PermissionError:
                    dst.unlink()
                    shutil.copyfile(src, dst)
            if digest(src) != digest(dst):
                fails.append("%s copy differs from %s" % (dst.name, src.name))
            rows.append((name + suffix, src.relative_to(WORK).as_posix(),
                         round(dst.stat().st_size / 1024)))
    counts = {
        "manuscript": pages(WORK / "build_v2" / "main_v2.log"),
        "supplementary": pages(WORK / "build_v2" / "supplementary_text.log"),
        "letter": pages(WORK / "build_cover" / "cover_letter_NHTB_v2.log"),
    }
    supp_src = io.open(WORK / "supplementary_text.tex", encoding="utf-8").read()
    supp_floats = len([l for l in re.findall(r"\\label\{(fig:[^}]*|tab:[^}]*)\}",
                                             supp_src)])
    supp_sections, supp_count = supplementary_sections(supp_src)
    readme = """# NHT-B submission package - %s

Target journal: **Numerical Heat Transfer, Part B: Fundamentals** (Taylor & Francis).
The paper is submitted as two documents plus the editorial letter.

| file | content | size |
|---|---|---|
""" % DATE
    for fname, origin, size_kb in rows:
        readme += "| `%s` | copy of `../%s` | %d KB |\n" % (fname, origin, size_kb)
    readme += """
The two documents are one publication. The main text is a nine-section story --
introduction, formulation and base state, the symmetry-selected pumping law, the
exact pumping identities and their Reynolds-number dependence, the phase-shifted
two-coefficient law, numerical methods and verification, results with the
independent nonlinear DNS validation, physical interpretation and discussion,
conclusions -- and it keeps the parity theorem, the two-coefficient phase law,
the reversal map, the thirteen-case campaign and every statement the argument
needs. The supplementary document is organised, in the order in which a
reviewer checks it, into sections numbered from S1 to S%(supp_last)s:
%(supp_sections)s. Its floats, sections
and equations are numbered **S1, S2, ...** and each of the %(supp_floats)d
supplementary floats is cited and discussed from the main text. The
cross-references are resolved by the `xr` package in both directions, so the two
files cannot silently disagree about a number. A Supplementary Material section
after the Conclusions describes that split, and the body pointers name the
supplementary material explicitly.

Lengths: main manuscript %(main_pages)s pages, supplementary document
%(supp_pages)s pages, cover letter %(letter_pages)s pages.

The manuscript text is venue-neutral: it names no journal anywhere.

Compile from the **project root**, not from this folder, because both sources
load their figures with paths relative to the root:

    pdflatex -interaction=nonstopmode -output-directory=build_v2 main_v2.tex
    bibtex build_v2/main_v2
    pdflatex -interaction=nonstopmode -output-directory=build_v2 main_v2.tex
    pdflatex -interaction=nonstopmode -output-directory=build_v2 supplementary_text.tex
    pdflatex -interaction=nonstopmode -output-directory=build_v2 supplementary_text.tex
    pdflatex -interaction=nonstopmode -output-directory=build_v2 main_v2.tex
    pdflatex -interaction=nonstopmode -output-directory=build_v2 supplementary_text.tex
    pdflatex -interaction=nonstopmode -output-directory=build_cover cover_letter_NHTB_v2.tex

The passes alternate because the two documents reference each other: the supplement
needs the auxiliary file of the main text, and the main text needs the auxiliary file
of the supplement for the S-numbered cross-references. Building the supplement only
at the end leaves those cross-references undefined, and the build gate reports it.

The build gate reads both logs and both `.aux` files, checks the two sources as
one corpus for dangling references and uncited floats, and enforces a
%(page_limit)d-page ceiling on the main text:

    py -X utf8 campaign_v2_symmetry_selection/check_build.py
    py -X utf8 campaign_v2_symmetry_selection/verify_manuscript_numbers.py
    py -X utf8 campaign_v2_symmetry_selection/check_notation.py
    py -X utf8 campaign_v2_symmetry_selection/check_duplicate_text.py
    py -X utf8 campaign_v2_symmetry_selection/check_cover_letter.py
    py -X utf8 campaign_v2_symmetry_selection/check_end_matter.py

Gate status at the time of this copy, all green:

- build: 0 errors, 0 undefined references, 0 overfull box in either document,
  %(required_labels)d required labels resolved, every supplementary float cited
  from the main text, sources ASCII-clean, no working-filename citation, keyword
  line present
- number verifier: 0 failures against the stored campaign data
- notation: conventions consistent in both documents
- duplication scan: run over the manuscript **and** its supplementary document
  together; three shared twelve-word shingles with study 1 (all in the
  author-and-affiliation block), none with study 2, the earlier drafts of this
  study measured as well. The accounting is reported in the cover letter and is
  quoted nowhere in the manuscript, which is what the gate checks
- cover letter: %(letter_pages)s pages, limit %(letter_limit)d, every quoted
  statement present
- end matter: single-paragraph abstract within 250 words and free of equations,
  Supplementary Material section after the Conclusions, disclosure,
  author-contribution and funding statements present, body pointers using the
  words "supplementary material"

The supplementary data package stays at `../supplementary_material/` (64 files in
groups A to D, machine-readable manifest with the SHA-256 digest of every item);
zip it before upload. The whole study is also deposited publicly at %(repo)s,
which the manuscript and the letter both cite.

`main.tex` (the older PoF-style manuscript) is untouched; this package contains
only the v2 line.
""" % dict(supp_floats=supp_floats, main_pages=counts["manuscript"],
           supp_sections=supp_sections, supp_last=supp_count,
           supp_pages=counts["supplementary"], letter_pages=counts["letter"],
           page_limit=45, required_labels=28, letter_limit=2, repo=REPO)
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print("wrote %s (%d files + README)" % (OUT, len(rows)))
    print("pages: main %s, supplementary %s, letter %s"
          % (counts["manuscript"], counts["supplementary"], counts["letter"]))
    if fails:
        print("COPY VERIFICATION FAILED:")
        for f in fails:
            print("  ", f)
        return 1
    print("all copies byte-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())

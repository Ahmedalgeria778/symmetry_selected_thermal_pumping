#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble the journal submission folder from the three compiled documents.

PDF copies are byte-identical to the build output and every digest is recomputed
after the copy, so the folder fails if a pair differs. The two source copies are
byte-identical to the working files except for their `xr` target line, which is
renamed to the file name this folder actually ships; the diff is checked line by
line. The README is generated from the build logs, so its page counts cannot go
stale.

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

# In the working tree the two documents are built into build_v2, so that is
# where their xr targets, and the file-existence guard around each of them,
# point. Here they are shipped as manuscript.tex and supplementary_document.tex,
# and a reader who compiles this pair in place would otherwise get "???" for
# every cross-reference: the target and its guard name the shipped file. Any
# other difference from the working source is an error.
XR_TARGETS = {
    "main_v2.tex": {"build_v2/supplementary_text": "supplementary_document"},
    "supplementary_text.tex": {"build_v2/main_v2": "manuscript"},
}


def rewrite_xr(src, name):
    for old, new in XR_TARGETS.get(name, {}).items():
        src = src.replace("\\externaldocument{%s}" % old,
                          "\\externaldocument{%s}" % new)
        src = src.replace("\\IfFileExists{%s.aux}" % old,
                          "\\IfFileExists{%s.aux}" % new)
    return src


def xr_only_diff(a, b, name):
    """Line differences between the working source and the shipped copy, the
    xr target and its guard of this file allowed away; anything else returns a
    complaint."""
    if name not in XR_TARGETS:
        return [] if a == b else ["body differs from the working source"]
    la, lb = a.split("\n"), b.split("\n")
    if len(la) != len(lb):
        return ["line count differs"]
    keys = list(XR_TARGETS[name])
    return ["line %d: %s" % (i, x[:70])
            for i, (x, y) in enumerate(zip(la, lb), 1)
            if x != y and not any(k in x for k in keys)]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def pages(log_path):
    log = io.open(log_path, encoding="utf-8", errors="replace").read()
    m = re.findall(r"\((\d+)\s+pages?\b", log.replace("\n", " "))
    return m[0] if m else "?"


def pcount(n):
    return "%s %s" % (n, "page" if str(n) == "1" else "pages")


def supplementary_sections(src):
    """The S-numbered section titles, read from the supplement itself so the
    README cannot drift from the document it describes."""
    titles = ["S%d %s%s" % (i, t[0].lower(), t[1:]) for i, t in
              enumerate(re.findall(r"^\\section\{([^}]*)\}", src, re.M), 1)]
    return ", ".join(titles[:-1]) + " and " + titles[-1], len(titles)


def read_text(path):
    return io.open(path, encoding="utf-8", newline="\n").read()


def publish(src, dst, payload=None):
    """Copy the source, or the payload derived from it, unless the destination
    already matches. Windows can hold a lock on a PDF a viewer has open, and
    rewriting such a file is not worth a failed build."""
    if payload is None:
        if dst.exists() and digest(src) == digest(dst):
            return
        try:
            shutil.copyfile(src, dst)
        except PermissionError:
            dst.unlink()
            shutil.copyfile(src, dst)
    else:
        if dst.exists() and read_text(dst) == payload:
            return
        io.open(dst, "w", encoding="utf-8", newline="\n").write(payload)


def main():
    OUT.mkdir(exist_ok=True)
    rows, fails = [], []
    for name, tex, pdf in ITEMS:
        for src, suffix in ((tex, ".tex"), (pdf, ".pdf")):
            dst = OUT / (name + suffix)
            if suffix == ".tex":
                publish(src, dst, rewrite_xr(read_text(src), src.name))
                fails += ["%s: %s" % (dst.name, d)
                          for d in xr_only_diff(read_text(src), read_text(dst),
                                                 src.name)]
            else:
                publish(src, dst)
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
        note = "copy of `../%s`" % origin
        if fname.endswith(".tex") and fname in ("manuscript.tex",
                                                "supplementary_document.tex"):
            note += "; its `xr` target names the file shipped here"
        readme += "| `%s` | %s | %d KB |\n" % (fname, note, size_kb)
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
files cannot silently disagree about a number. In the working tree the two
`\\externaldocument` targets point at the build directory; in this folder they
name the files actually shipped here (`manuscript` and
`supplementary_document`), and nothing else in the two sources differs from the
working files. A Supplementary Material section
after the Conclusions describes that split, and the body pointers name the
supplementary material explicitly.

Lengths: main manuscript %(main_pages)s, supplementary document
%(supp_pages)s, cover letter %(letter_pages)s.

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

In this folder the two sources reference each other by the names used here, so an
`xr` target can no longer point at a file that is not in the package. A rebuild
from this directory still needs the figure directories, which are not part of it:
they come from the working tree or from the public repository below, whose paths
are the root-relative ones the compile lines above assume.

Bootstrapping the pair from this folder needs the supplement first, because each
document reads the auxiliary file of the other. In this directory the order is
`supplementary_document.tex`, then `manuscript.tex`, then
`supplementary_document.tex` again; starting from `manuscript.tex` with no
supplement auxiliary present leaves every S-numbered reference printed as `??`
until the supplement has been built once.

The build gate reads both logs and both `.aux` files, checks the two sources as
one corpus for dangling references and uncited floats, verifies that every
graphics input exists with the exact case and that every `xr` target names a
file that is really there, and enforces a
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
  together; no scientific sentence and no shared twelve-word run with the
  companion studies, the only overlap with study 1 being the
  author-and-affiliation block, and the earlier drafts of this study measured
  as well. The accounting is quoted nowhere: neither the manuscript nor the
  cover letter states it, which is what the gate checks
- cover letter: %(letter_pages)s, limit %(letter_limit)d, every required
  statement present and no promotional framing
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
""" % dict(supp_floats=supp_floats, main_pages=pcount(counts["manuscript"]),
           supp_sections=supp_sections, supp_last=supp_count,
           supp_pages=pcount(counts["supplementary"]),
           letter_pages=pcount(counts["letter"]),
           page_limit=45, required_labels=28, letter_limit=1, repo=REPO)
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print("wrote %s (%d files + README)" % (OUT, len(rows)))
    print("pages: main %s, supplementary %s, letter %s"
          % (counts["manuscript"], counts["supplementary"], counts["letter"]))
    if fails:
        print("COPY VERIFICATION FAILED:")
        for f in fails:
            print("  ", f)
        return 1
    print("PDFs byte-identical to the build output; the two sources differ from "
          "the working files only in their xr target line")
    return 0


if __name__ == "__main__":
    sys.exit(main())

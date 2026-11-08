"""Plan item F.2: word-level duplication check of main_v2.tex against the
companion manuscripts of the mixed-convection series.

The comparison is deliberately strict: LaTeX is stripped to plain words, then
shared sentences (exact) and shared word-shingles (near-exact) are reported so
that any carried-over prose can be found and rewritten before submission.
The titles and abstracts of every companion are printed first, because the
series membership of each file has to be established rather than assumed.
"""
import glob
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHINGLE = 12


def strip_latex(text):
    text = re.sub(r"(?<!\\)%.*", "", text)
    text = re.sub(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", r" \1 ", text,
                  flags=re.S)
    text = re.sub(r"\\\\\[(.*?)\\\\\]", " ", text, flags=re.S)
    text = re.sub(r"\$(.*?)\$", " ", text, flags=re.S)
    text = re.sub(r"\\[a-zA-Z]+\*?", " ", text)
    text = re.sub(r"[{}()\[\]&]", " ", text)
    return re.sub(r"\s+", " ", text)


def words(text):
    return re.findall(r"[a-z]{3,}", text.lower())


def sentences(text):
    return [s.strip() for s in re.split(r"(?<=[.;:]) ", text) if len(s) > 60]


def title_and_abstract(fname):
    raw = io.open(fname, encoding="utf-8", errors="replace").read()
    t = re.search(r"\\title\{(.{0,400}?)\}", raw, re.S)
    a = re.search(r"\\begin\{abstract\}(.{0,2000}?)\\end\{abstract\}", raw, re.S)
    tt = re.sub(r"\s+", " ", t.group(1)) if t else "(no \\title)"
    aa = re.sub(r"\s+", " ", a.group(1))[:400] if a else "(no abstract)"
    return tt, aa


def pdf_words(path):
    from pypdf import PdfReader
    txt = " ".join((p.extract_text() or "") for p in PdfReader(path).pages)
    return words(re.sub(r"\s+", " ", txt))


def pdf_title(path):
    """First 110 characters of page one, to prove which version was scanned."""
    from pypdf import PdfReader
    txt = re.sub(r"\s+", " ", (PdfReader(path).pages[0].extract_text() or ""))
    return txt[:110]


def longest_runs(mw, other_shingles, k=4):
    """Longest contiguous stretches of the manuscript whose word shingles also
    occur in the other document."""
    runs = []
    i = 0
    while i + SHINGLE <= len(mw):
        if tuple(mw[i:i + SHINGLE]) in other_shingles:
            j = i
            while (j + SHINGLE <= len(mw)
                   and tuple(mw[j:j + SHINGLE]) in other_shingles):
                j += 1
            runs.append((i, j + SHINGLE - i))
            i = j
        else:
            i += 1
    runs.sort(key=lambda r: -r[1])
    return runs[:k]


def main():
    main_tex = os.path.join(ROOT, "main_v2.tex")
    supp_tex = os.path.join(ROOT, "supplementary_text.tex")
    companions = sorted(
        p for p in glob.glob(os.path.join(ROOT, "*.tex"))
        if os.path.basename(p) not in ("main_v2.tex", "supplementary_text.tex",
                                       "main.tex",
                                       "table2_corrected.tex",
                                       "table_PiT_k10_vs_Pe.tex"))
    # companion papers of the series that exist only as PDFs. Study 1 is checked
    # against BOTH the project copy and the final version kept outside the tree,
    # so the footnote cannot be validated against a stale draft.
    pdfs = {
        "paper1_nonmonotonic_integral (NHT-B)":
            os.path.join(os.path.dirname(ROOT), "Paper_developpement",
                         "Paper_Final_07_08_2026.pdf"),
        "paper1_nonmonotonic_integral (final version, Downloads)":
            os.path.join(os.path.expanduser("~"), "Downloads",
                         "Article_BENIAICHE.pdf"),
        "paper2_semi_analytical (PoF)":
            os.path.join(ROOT, "Semi_analytical_Solution_laminar_mixed_convection.pdf"),
        "neural_operator_companion":
            os.path.join(os.path.dirname(ROOT), "Paper_developpement",
                         "Main_manuscript_revised.pdf"),
    }
    study1 = [k for k in pdfs if k.startswith("paper1_nonmonotonic")]
    study2 = [k for k in pdfs if k.startswith("paper2_semi")]
    print("== series candidates ==")
    for c in companions:
        tt, aa = title_and_abstract(c)
        print("\n%s\n  title: %s\n  abstract: %s" % (os.path.basename(c), tt, aa))

    # The publication is the two-document set: the scan covers the manuscript
    # and its supplementary document together, so a shared phrase cannot hide
    # in the part that is not scanned.
    corpus = (io.open(main_tex, encoding="utf-8", errors="replace").read() + "\n"
              + io.open(supp_tex, encoding="utf-8", errors="replace").read())
    mw = words(strip_latex(corpus))
    ms = set(" ".join(s.split()) for s in sentences(strip_latex(corpus)))
    m_sh = set(tuple(mw[i:i + SHINGLE]) for i in range(len(mw) - SHINGLE))
    print("\nmanuscript plus supplementary document: %d words, %d sentences, "
          "%d %d-word shingles"
          % (len(mw), len(ms), len(m_sh), SHINGLE))

    worst = []
    res = {}
    for label, path in sorted(pdfs.items()):
        if not os.path.exists(path):
            print("\n== %s : NOT FOUND (%s) ==" % (label, path))
            continue
        cw = pdf_words(path)
        c_sh = set(tuple(cw[i:i + SHINGLE]) for i in range(len(cw) - SHINGLE))
        shared_sh = m_sh & c_sh
        frac = (len(shared_sh) / len(m_sh)) if m_sh else 0.0
        print("\n== %s ==" % label)
        print("  scanned document : %s" % pdf_title(path))
        print("  shared %d-word shingles : %d (%.3f%% of manuscript)"
              % (SHINGLE, len(shared_sh), 100 * frac))
        res[label] = len(shared_sh)
        for start, ln in longest_runs(mw, c_sh):
            print("    RUN(%d words): %s" % (ln, " ".join(mw[start:start + ln])[:170]))
            worst.append((label, ln, " ".join(mw[start:start + ln])))

    for c in companions:
        txt = strip_latex(io.open(c, encoding="utf-8", errors="replace").read())
        cw = words(txt)
        cs = set(" ".join(s.split()) for s in sentences(txt))
        c_sh = set(tuple(cw[i:i + SHINGLE]) for i in range(len(cw) - SHINGLE))
        shared_sent = ms & cs
        shared_sh = m_sh & c_sh
        frac = (len(shared_sh) / len(m_sh)) if m_sh else 0.0
        print("\n== %s ==" % os.path.basename(c))
        print("  exact shared sentences : %d" % len(shared_sent))
        print("  shared %d-word shingles  : %d (%.3f%% of manuscript)"
              % (SHINGLE, len(shared_sh), 100 * frac))
        res[os.path.basename(c)] = 100 * frac
        for s in sorted(shared_sent)[:6]:
            print("    SENT: %s" % s[:150])
        runs = longest_runs(mw, c_sh, k=5)
        for start, ln in runs:
            print("    RUN(%d words): %s" % (ln, " ".join(mw[start:start + ln])[:170]))
            worst.append((os.path.basename(c), ln, " ".join(mw[start:start + ln])))

    print("\n== longest duplicated runs across all companions ==")
    for name, ln, frag in sorted(worst, key=lambda r: -r[1])[:8]:
        print("  %-38s %3d words: %s" % (name, ln, frag[:150]))

    # The review asked that the duplication accounting live in the cover
    # letter, not in a manuscript footnote: the letter must quote what this
    # scan measures, and the manuscript must quote nothing of the sort.
    flat = re.sub(r"\s+", " ", io.open(main_tex, encoding="utf-8",
                                       errors="replace").read())
    letter = re.sub(r"\s+", " ", io.open(
        os.path.join(os.path.dirname(main_tex), "cover_letter_NHTB_v2.tex"),
        encoding="utf-8", errors="replace").read())
    fails = []
    m = re.search(r"returns ([a-z0-9]+) shared twelve-word shingles", letter)
    first = max([res[k] for k in study1 if k in res], default=None)
    second = max([res[k] for k in study2 if k in res], default=None)
    if len([k for k in study1 if k in res]) < len(study1):
        fails.append("not every study-1 version was scanned: %s"
                     % [k for k in study1 if k not in res])
    spelled = {"none": 0, "zero": 0, "one": 1, "two": 2, "three": 3, "four": 4}
    quoted_first = spelled.get(m.group(1).lower(), None) if m else None
    if quoted_first is None and m is not None:
        quoted_first = int(m.group(1))
    if m is None:
        fails.append("duplication sentence not found in the cover letter")
    elif quoted_first != first or second != 0:
        fails.append("cover letter quotes %s shingles, measured %s against study 1 "
                     "and %s against the unmentioned second companion, which must "
                     "stay at zero" % (m.group(1), first, second))
    else:
        print("\ncover-letter companion-study shingles : %d  ok  "
              "(study 1 scanned in %d version(s); the unsubmitted companion measures %d "
              "and is not named in the text)"
              % (first, len([k for k in study1 if k in res]), second))
    for name, doc in (("manuscript", flat), ("cover letter", letter)):
        if re.search(r"twelve-word shingle|draft overlap|earlier draft", doc):
            if name == "manuscript":
                fails.append("the manuscript still carries duplication-accounting "
                             "prose, which the review asked to delete")
    drafts = sorted((v, k) for k, v in res.items()
                    if k.endswith(".tex") and k != "cover_letter_NHTB_v2.tex")
    if len(drafts) >= 2:
        top = [round(v, 1) for v, _ in sorted(drafts, reverse=True)[:2]]
        print("draft overlaps, measured only (quoted nowhere): %s (%s)"
              % (top, ",".join(k for _, k in sorted(drafts, reverse=True)[:2])))
    for f in fails:
        print("DUPLICATION MISMATCH: " + f)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

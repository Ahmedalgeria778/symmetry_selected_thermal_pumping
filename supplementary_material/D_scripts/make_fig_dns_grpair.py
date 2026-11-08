"""Figure for the nonlinear (DNS) test of Corollary 1.

Corollary 1 states Pi_T = Gr * A(k, Re, Pr), so the reduced coefficient A must be
a single curve for every Grashof number, and its zero must be the reversal point.
The figure draws the two perturbative branches at Gr = 50 and Gr = 200 and
overlays the four DNS measurements of the Grashof-pair matrix, so that the
nonlinear confirmation is read as points falling on one line rather than as a
table of ratios.

    py -X utf8 campaign_v2_symmetry_selection/make_fig_dns_grpair.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
WORK = OUT.parent
sys.path.insert(0, str(WORK))
_spec = importlib.util.spec_from_file_location(
    "d345", OUT / "run_dns_D3_D4_D5.py")
d345 = importlib.util.module_from_spec(_spec)
sys.modules["d345"] = d345
_spec.loader.exec_module(d345)

K, PE = 0.8, 7.0
GR_LO, GR_HI = 50.0, 200.0
MEASURED = OUT / "DNS_D2_grpair_measured.csv"
RE_C = OUT / "D2_Re_c_vs_Gr.csv"
PNG = OUT / "FIG_v2_11_dns_grpair.png"
JSON = OUT / "FIG_v2_11_dns_grpair.json"


def main():
    df = pd.read_csv(MEASURED)
    rec = pd.read_csv(RE_C)
    re_c = float(rec.Re_c.iloc[0])

    grid = np.linspace(12.0, 60.0, 25)
    curves = {}
    for Gr in (GR_LO, GR_HI):
        curves[Gr] = np.array([d345.coefficients(K, Re, PE, Gr)[0] / Gr
                               for Re in grid])
        print("  perturbative branch Gr=%5.1f computed at %d Reynolds numbers"
              % (Gr, len(grid)), flush=True)
    collapse = float(np.max(np.abs(curves[GR_LO] - curves[GR_HI])
                            / np.abs(curves[GR_LO])))
    print("  max relative separation of the two branches: %.3e" % collapse,
          flush=True)

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.plot(grid, curves[GR_LO], "-", color="tab:blue", lw=1.4,
            label=r"perturbative, $\mathrm{Gr}=50$")
    ax.plot(grid, curves[GR_HI], "--", color="tab:red", lw=1.4,
            label=r"perturbative, $\mathrm{Gr}=200$")
    marks = {GR_LO: ("o", "tab:blue", -0.7), GR_HI: ("s", "tab:red", 0.7)}
    for row in df.itertuples():
        m, col, dx = marks[float(row.Gr)]
        ax.plot(row.Re + dx, row.A_dns, m, color=col, ms=8, mec="k", mew=0.7,
                zorder=5)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axvline(re_c, color="gray", lw=0.8, ls=":")
    ax.text(re_c, ax.get_ylim()[1], r" $Re_c=%.2f$" % re_c, ha="left",
            va="top", fontsize=8, color="gray")
    ax.set_xlabel(r"$Re$")
    ax.set_ylabel(r"$A=\Pi_T/\mathrm{Gr}$   ($k=0.8$, $Pe=7$)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("DNS Grashof pair on the single perturbative branch "
                 "(Corollary 1)", fontsize=10)
    fig.tight_layout()
    fig.savefig(PNG, dpi=200)
    plt.close(fig)

    agree = {}
    for Re in sorted(set(df.Re)):
        a = df[np.isclose(df.Re, Re)].sort_values("Gr").A_dns.values
        agree[str(int(Re))] = float(abs(a[1] - a[0]) / abs(a.mean()))
    summary = dict(k=K, Pe=PE, Re_c=re_c, branch_collapse_rel=collapse,
                   marker_offset_re=0.7, pair_agreement_rel=agree,
                   grid_points=int(len(grid)),
                   dns_points=[dict(Re=float(r.Re), Gr=float(r.Gr),
                                    A_dns=float(r.A_dns),
                                    A_perturbative=float(r.A_perturbative),
                                    A_rel_err=float(r.A_rel_err))
                               for r in df.itertuples()])
    JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("saved %s and %s" % (PNG.name, JSON.name))
    for Re, v in agree.items():
        print("  pair agreement in A at Re=%s: %.3e relative" % (Re, v))
    for p in summary["dns_points"]:
        print("  Re=%5.1f Gr=%5.1f  A_dns=%.6e  A_per=%.6e  rel=%.3e"
              % (p["Re"], p["Gr"], p["A_dns"], p["A_perturbative"],
                 p["A_rel_err"]))


if __name__ == "__main__":
    main()

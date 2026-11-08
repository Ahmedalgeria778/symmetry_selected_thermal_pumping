#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extend the D.1 Re_c(k) cartography to the two long-wavelength numbers.

The D.1 scan stopped at Re = 40, so k = 0.4 and k = 0.6 were reported as
"single sign over the scanned interval" (>$40$ in Table tab:Rec_k). Brent
refinement on a widened interval shows both wavenumbers do reverse, at
Re = 79.08 and 51.25: the reversal locus is defined at all seven wavenumbers and
the empirical power law can be re-fitted on seven points instead of five.

Outputs
    D1_Re_grid_lowk.csv     extended grid (same columns as D1_Re_grid.csv)
    D1_Re_c_all.csv         merged seven-point refined locus
    D1_Re_c_fit.json        power-law fit and per-point deviations
    FIG_v2_02_Re_c_vs_k.png regenerated from the merged locus
    FIG_v2_01_piT_vs_Re_kartsch.png regenerated from the merged grid
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

WORK = Path(r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification")
OUT = WORK / "campaign_v2_symmetry_selection"
_spec = importlib.util.spec_from_file_location(
    "v4", WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py")
v4 = importlib.util.module_from_spec(_spec)
sys.modules["v4"] = v4
_spec.loader.exec_module(v4)

PE, GR = 7.0, 100.0
RES = (5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
       60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0)
KS = (0.4, 0.6)


def pi_c(k, Re):
    """Second-order mixed pumping coefficient at unit symmetric amplitude."""
    ts = v4.solve_thermal_bvp(k, PE, True)
    ta = v4.solve_thermal_bvp(k, PE, False)
    hs = v4.solve_hydro_bvp(k, Re, GR, ts)
    ha = v4.solve_hydro_bvp(k, Re, GR, ta)
    y = np.linspace(-1.0, 1.0, 1001)
    Ts, Ta = ts.theta(y), ta.theta(y)
    Ts_y, Ta_y = ts.theta_y(y), ta.theta_y(y)
    Us, Ua, Vs, Va = hs.U(y), ha.U(y), hs.V(y), ha.V(y)
    Z1 = Ua * np.conjugate(-1j * k * Ts) + Va * np.conjugate(Ts_y)
    Z2 = Us * np.conjugate(-1j * k * Ta) + Vs * np.conjugate(Ta_y)
    Kc = 0.5 * (np.real(Z1) + np.real(Z2))
    m = v4.solve_mean_bvp(y, Kc, PE)
    thy = m.theta_y(y)
    return 0.5 * (-thy[-1] - thy[0]), abs(-thy[-1] + thy[0])


def main():
    rows = []
    for k in KS:
        for Re in RES:
            if Re <= 40.0 and k in (0.4, 0.6):
                continue          # already in D1_Re_grid.csv, do not duplicate
            val, anti = pi_c(k, Re)
            rows.append(dict(k=k, Pe=PE, Gr=GR, Re=Re, Pr=PE / Re,
                             Pi_T=val, A=val / GR, antisym_err=anti))
            print("k=%.1f Re=%5.1f  Pi_T=%+.8e  (asym err %.1e)"
                  % (k, Re, val, anti), flush=True)
    gdf = pd.DataFrame(rows)
    gdf.to_csv(OUT / "D1_Re_grid_lowk.csv", index=False)

    merged_grid = pd.concat(
        [pd.read_csv(OUT / "D1_Re_grid.csv"), gdf], ignore_index=True)

    loc = []
    for k in (0.4, 0.6):
        sub = merged_grid[merged_grid.k == k].sort_values("Re")
        lo = sub[sub.Pi_T < 0].Re.max()
        hi = sub[sub.Pi_T > 0].Re.min()
        f = lambda R, k=k: pi_c(k, R)[0]
        rec_c = brentq(f, lo, hi, xtol=1e-8, rtol=1e-10, maxiter=200)
        loc.append(dict(k=k, Re_c=rec_c, bracket_lo=lo, bracket_hi=hi,
                        Pi_at_lo=f(lo), Pi_at_hi=f(hi)))
    for r in loc:
        print("Re_c(k=%.2f) = %.9f from bracket (%.1f, %.1f)"
              % (r["k"], r["Re_c"], r["bracket_lo"], r["bracket_hi"]))

    old = pd.read_csv(OUT / "D1_Re_c.csv")
    cdf = pd.concat([old[~old.k.isin(KS)], pd.DataFrame(loc)],
                    ignore_index=True).sort_values("k")
    cdf.to_csv(OUT / "D1_Re_c_all.csv", index=False)

    lk, lr = np.log(cdf.k.values), np.log(cdf.Re_c.values)
    n, lnC = np.polyfit(lk, lr, 1)
    C = math.exp(lnC)
    pred = C * cdf.k.values ** n
    dev = np.abs(pred - cdf.Re_c.values) / cdf.Re_c.values
    five = cdf[cdf.k >= 0.8]
    n5, l5 = np.polyfit(np.log(five.k), np.log(five.Re_c), 1)
    C5 = math.exp(l5)
    p5 = C5 * five.k.values ** n5
    dev5 = np.abs(p5 - five.Re_c.values) / five.Re_c.values

    fit = dict(form="Re_c = C k^n", C=C, n=n, n_points=len(cdf),
               k_list=[float(v) for v in cdf.k],
               Re_c_list=[float(v) for v in cdf.Re_c],
               max_rel_deviation=float(dev.max()),
               rel_deviations=[float(v) for v in dev],
               five_point=dict(C=C5, n=n5, n_points=len(five),
                               k_list=[float(v) for v in five.k],
                               max_rel_deviation=float(dev5.max())))
    json.dump(fit, open(OUT / "D1_Re_c_fit.json", "w"), indent=2)
    print("\nseven-point fit: Re_c = %.4f k^(%.4f), max rel dev %.3f%%"
          % (C, n, 100 * dev.max()))
    for k, rc, p, d in zip(cdf.k, cdf.Re_c, pred, dev):
        print("   k=%.1f  Re_c=%9.4f  fit=%9.4f  dev=%.2f%%" % (k, rc, p, 100 * d))
    print("five-point fit (k >= 0.8): Re_c = %.4f k^(%.4f), max rel dev %.3f%%"
          % (C5, n5, 100 * dev5.max()))

    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.plot(cdf.k, cdf.Re_c, "s-", color="tab:red")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$Re_c$  ($\Pi_T=0$, $Pe=7$, $Gr=100$)")
    for _, r in cdf.iterrows():
        ax.annotate("%.2f" % r.Re_c, (r.k, r.Re_c), textcoords="offset points",
                    xytext=(6, 4), fontsize=8)
    kk = np.linspace(cdf.k.min(), cdf.k.max(), 200)
    ax.plot(kk, C * kk ** n, "--", color="tab:blue", lw=1.0,
            label=r"fit $\mathrm{Re_c}=%.2f\,k^{%.2f}$" % (C, n))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_02_Re_c_vs_k.png", dpi=200)
    print("\nsaved FIG_v2_02_Re_c_vs_k.png (seven points)")

    # FIG_v2_01 regenerated on the merged grid so that every curve contains
    # its own crossing; the k >= 0.8 curves still stop at Re = 40, which is
    # where their scanned interval ends.
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for k in sorted(merged_grid.k.unique()):
        sub = merged_grid[merged_grid.k == k].sort_values("Re")
        ax.plot(sub.Re, sub.A, "o-", ms=4, lw=1.2, label="k = %.1f" % k)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"$Re$ at fixed $Pe=7$")
    ax.set_ylabel(r"$A=\Pi_T/\mathrm{Gr}$")
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.legend(fontsize=8, ncol=2)
    ax.set_title(r"Symmetry-selected pumping: reversal vs $Re$ for seven $k$ (perturbative)")
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_01_piT_vs_Re_kartsch.png", dpi=200)
    plt.close(fig)
    print("saved FIG_v2_01_piT_vs_Re_kartsch.png (merged grid)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

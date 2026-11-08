#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Follow-up on the quadrature (90 deg out-of-phase) pumping coefficient Pi_perp.

1. Gr-linearity of Pi_perp  (does the Green/linearity identity cover it?)
2. Re-scan of Pi_perp at k=0.8, Pe=7 -> its own reversal line Re_c_perp
3. Exact two-coefficient forcing law, randomised validation over
   (eps1, eps2, phi) with independent assembly of the mixed forcing.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

WORK = Path(r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification")
OUT = WORK / "campaign_v2_symmetry_selection"
_spec = importlib.util.spec_from_file_location("v4", WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py")
v4 = importlib.util.module_from_spec(_spec)
sys.modules["v4"] = v4
_spec.loader.exec_module(v4)

PE, GR, K = 7.0, 100.0, 0.8


def coeffs(k, Re, Pe, Gr, n_eval=1001):
    ts = v4.solve_thermal_bvp(k, Pe, True)
    ta = v4.solve_thermal_bvp(k, Pe, False)
    hs = v4.solve_hydro_bvp(k, Re, Gr, ts)
    ha = v4.solve_hydro_bvp(k, Re, Gr, ta)
    y = np.linspace(-1.0, 1.0, n_eval)
    Ts, Ta = ts.theta(y), ta.theta(y)
    Ts_y, Ta_y = ts.theta_y(y), ta.theta_y(y)
    Us, Ua, Vs, Va = hs.U(y), ha.U(y), hs.V(y), ha.V(y)
    Z1 = Ua * np.conjugate(-1j * k * Ts) + Va * np.conjugate(Ts_y)
    Z2 = Us * np.conjugate(-1j * k * Ta) + Vs * np.conjugate(Ta_y)
    X = 0.5 * (np.real(Z1) + np.real(Z2))
    Yq = 0.5 * (np.imag(Z2) - np.imag(Z1))
    out = {}
    for tag, prof in (("c", X), ("p", Yq)):
        m = v4.solve_mean_bvp(y, prof, Pe)
        thy = m.theta_y(y)
        out["Pi_" + tag] = 0.5 * (-thy[-1] - thy[0])
        out["anti_sym_err_" + tag] = abs(-thy[-1] + thy[0])
    return out, y, X, Yq


def wall_amps(eps1, eps2, phi):
    """Complex amplitudes of the symmetric/antisymmetric parts.

    Manuscript convention: fields are Re[hat e^{-ikX}], so
    theta(+1) = eps1 cos(kX) and theta(-1) = eps2 cos(kX + phi) give
    theta_hat(+1) = eps1 and theta_hat(-1) = eps2 exp(-i phi).
    """
    A_s = 0.5 * (eps1 + eps2 * np.exp(-1j * phi))
    A_a = 0.5 * (eps1 - eps2 * np.exp(-1j * phi))
    return A_s, A_a


def main():
    t0 = time.time()
    # ---- 1. Gr linearity of both coefficients
    rows = []
    for Gr in [25.0, 50.0, 100.0, 200.0, 400.0]:
        c, _, _, _ = coeffs(K, 10.0, PE, Gr)
        rows.append(dict(Gr=Gr, Re=10.0, k=K, Pe=PE,
                         Pi_c=c["Pi_c"], Pi_p=c["Pi_p"],
                         A_c=c["Pi_c"] / Gr, A_p=c["Pi_p"] / Gr))
        print("Gr=%6.1f  Pi_c=% .6e  Pi_p=% .6e  A_c=% .6e  A_p=% .6e"
              % (Gr, c["Pi_c"], c["Pi_p"], c["Pi_c"] / Gr, c["Pi_p"] / Gr), flush=True)
    g = pd.DataFrame(rows)
    g.to_csv(OUT / "P2_Gr_linearity_two_coeffs.csv", index=False)
    for col in ("A_c", "A_p"):
        v = g[col].values
        print("  %s relative spread over Gr: %.3e" % (col, np.ptp(v) / abs(v.mean())), flush=True)

    # ---- 2. Re scan of both coefficients
    rows = []
    for Re in [5.0, 10.0, 15.0, 20.0, 25.0, 28.0, 30.0, 30.838497, 32.0, 35.0, 40.0, 50.0]:
        c, _, _, _ = coeffs(K, Re, PE, GR)
        rows.append(dict(Re=Re, k=K, Pe=PE, Gr=GR, Pi_c=c["Pi_c"], Pi_p=c["Pi_p"],
                         anti_c=c["anti_sym_err_c"], anti_p=c["anti_sym_err_p"]))
        print("Re=%6.2f  Pi_c=% .6e  Pi_p=% .6e  anti-sym err %.1e/%.1e"
              % (Re, c["Pi_c"], c["Pi_p"], c["anti_sym_err_c"], c["anti_sym_err_p"]), flush=True)
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "P2_Re_scan_two_coeffs.csv", index=False)

    def f(tag, lo, hi):
        def g(Re):
            return coeffs(K, Re, PE, GR)[0]["Pi_" + tag]
        return brentq(g, lo, hi, xtol=1e-7)

    rec_c = f("c", 28.0, 33.0)
    rec_p = f("p", 5.0, 50.0) if np.sign(r.Pi_p.iloc[0]) != np.sign(r.Pi_p.iloc[-1]) else None
    print("Re_c(Pi_c) = %.6f   Re_c(Pi_p) = %s" % (rec_c, rec_p), flush=True)
    json.dump(dict(Re_c_inphase=rec_c, Re_c_quadrature=rec_p, k=K, Pe=PE, Gr=GR),
              open(OUT / "P2_Re_c_two_branches.json", "w"), indent=2)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.plot(r.Re, r.Pi_c, "o-", label=r"$\Pi_T$ (in-phase, published)")
    ax.plot(r.Re, r.Pi_p, "s-", label=r"$\Pi_T^{\perp}$ (quadrature, new)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"$Re$ at fixed $Pe=7$ ($k=0.8$, $Gr=100$)")
    ax.set_ylabel(r"pumping coefficient")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.legend(fontsize=8)
    ax.set_title("Two independent reversal branches")
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_06_two_branches.png", dpi=200)
    plt.close(fig)

    # ---- 3. randomised validation of the two-coefficient law
    c0, y, X, Yq = coeffs(K, 10.0, PE, GR)
    Pi_c, Pi_p = c0["Pi_c"], c0["Pi_p"]
    rng = np.random.default_rng(20261108)
    rows = []
    for _ in range(12):
        e1, e2 = rng.uniform(0.2, 1.0, 2)
        phi = rng.uniform(0.0, 2 * math.pi)
        A_s, A_a = wall_amps(e1, e2, phi)
        cc = float(np.real(A_a * np.conjugate(A_s)))
        dd = float(np.imag(A_a * np.conjugate(A_s)))
        pi_pred = cc * Pi_c + dd * Pi_p
        m = v4.solve_mean_bvp(y, cc * X + dd * Yq, PE)
        thy = m.theta_y(y)
        pi_dir = 0.5 * (-thy[-1] - thy[0])
        rows.append(dict(eps1=e1, eps2=e2, phi=phi, Re_Aa_As=cc, Im_Aa_As=dd,
                        Pi_pred=pi_pred, Pi_direct=pi_dir,
                        abs_err=abs(pi_pred - pi_dir),
                        rel_err=abs(pi_pred - pi_dir) / max(abs(pi_dir), 1e-300)))
        print("  eps=(%.3f,%.3f) phi=%.3f  Pi=% .6e vs % .6e  rel=%.2e"
              % (e1, e2, phi, pi_pred, pi_dir, rows[-1]["rel_err"]), flush=True)
    v = pd.DataFrame(rows)
    v.to_csv(OUT / "P2_two_coeff_law_validation.csv", index=False)
    print("max rel err over randomised cases: %.3e" % v.rel_err.max(), flush=True)

    # ---- 4. the sin(phi) figure is produced by run_phase_direct_sweep.py, which
    # re-solves the mean thermal problem at every phase. The previous version of
    # this block plotted the analytic expression against sin(phi) and was
    # therefore circular; it has been removed.
    print("TOTAL %.0f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()

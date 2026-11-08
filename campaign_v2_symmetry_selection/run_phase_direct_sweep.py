#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Honest phase sweep for the sin(phi) selection rule.

The previous version of FIG_v2_07_sin_phi.png plotted the analytic expression
0.5*eps^2*sin(phi)*Pi_perp against sin(phi), i.e. the identity against itself.
Here every point is obtained by assembling the mixed forcing for the requested
phase relation, solving the second-order mean thermal BVP numerically and
reading the wall-flux redistribution off that solution. The comparison with
sin(phi) is then a real test.

Outputs: D5_phase_sweep_direct.csv and FIG_v2_07_sin_phi.png
"""
from __future__ import annotations

import importlib.util
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = Path(r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification")
OUT = WORK / "campaign_v2_symmetry_selection"
_spec = importlib.util.spec_from_file_location(
    "v4", WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py")
v4 = importlib.util.module_from_spec(_spec)
sys.modules["v4"] = v4
_spec.loader.exec_module(v4)

K, PE, GR, RE = 0.8, 7.0, 100.0, 10.0
EPS = 1.0
N_PHI = 25


def kernels(k, Re, Pe, Gr, n_eval=1001):
    """Unit-amplitude kernels K_c, K_q and the corresponding coefficients."""
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
    Kc = 0.5 * (np.real(Z1) + np.real(Z2))
    Kq = 0.5 * (np.imag(Z2) - np.imag(Z1))
    return y, Kc, Kq


def wall_flux(y, profile, Pe):
    m = v4.solve_mean_bvp(y, profile, Pe)
    thy = m.theta_y(y)
    return 0.5 * (-thy[-1] - thy[0]), abs(-thy[-1] + thy[0])


def main():
    t0 = time.time()
    y, Kc, Kq = kernels(K, RE, PE, GR)
    Pi_c, anti_c = wall_flux(y, Kc, PE)
    Pi_p, anti_p = wall_flux(y, Kq, PE)
    print("Pi_c = % .8e (antisym err %.1e)" % (Pi_c, anti_c))
    print("Pi_p = % .8e (antisym err %.1e)" % (Pi_p, anti_p), flush=True)

    phis = np.linspace(0.0, math.pi, N_PHI)
    rows = []
    for phi in phis:
        # Manuscript convention: fields are Re[hat e^{-ikX}], so the wall data are
        # theta_hat(+1) = eps1, theta_hat(-1) = eps2 e^{-i phi}, i.e.
        # theta(X,+1) = eps1 cos(kX), theta(X,-1) = eps2 cos(kX + phi).
        A_s = 0.5 * EPS * (1.0 + np.exp(-1j * phi))
        A_a = 0.5 * EPS * (1.0 - np.exp(-1j * phi))
        c = float(np.real(A_a * np.conjugate(A_s)))
        d = float(np.imag(A_a * np.conjugate(A_s)))
        prof = c * Kc + d * Kq
        pi_dir, anti = wall_flux(y, prof, PE)
        pi_pred = c * Pi_c + d * Pi_p
        rows.append(dict(phi=phi, phi_over_pi=phi / math.pi, eps=EPS, k=K, Re=RE,
                         Pe=PE, Gr=GR, Re_Aa_As=c, Im_Aa_As=d,
                         Nu_direct=pi_dir, Nu_algebraic=pi_pred,
                         antisym_err=anti, sin_phi=math.sin(phi)))
        print("phi/pi=%.4f  c=%+.3e d=%+.3e  Nu=% .6e  pred=% .6e  err=%.1e"
              % (phi / math.pi, c, d, pi_dir, pi_pred, abs(pi_dir - pi_pred)),
              flush=True)

    df = pd.DataFrame(rows)
    ref = df.loc[(df.phi_over_pi - 0.5).abs().argsort().iloc[0], "Nu_direct"]
    df["Nu_ratio"] = df.Nu_direct / ref
    df["dev_from_sin"] = df.Nu_ratio - df.sin_phi
    df["abs_err_pred"] = (df.Nu_direct - df.Nu_algebraic).abs()
    df.to_csv(OUT / "D5_phase_sweep_direct.csv", index=False)

    print("\nmax |Nu(phi)/Nu(pi/2) - sin(phi)| = %.3e" % df.dev_from_sin.abs().max())
    print("max |Nu_direct - Nu_algebraic|      = %.3e" % df.abs_err_pred.max())
    print("max antisymmetry residual           = %.3e" % df.antisym_err.max())
    print("c vanishes to                        = %.3e" % df.Re_Aa_As.abs().max())

    fig, ax = plt.subplots(figsize=(5.8, 4.2), dpi=200)
    ax.plot(df.phi_over_pi, df.Nu_ratio, "o-", color="#1f77b4", ms=5,
            label="computed wall flux (mean BVP solved per phase)")
    ax.plot(df.phi_over_pi, df.sin_phi, "k--", lw=1.4,
            label=r"$\sin\phi$")
    ax.set_xlabel(r"$\phi/\pi$   (equal wall amplitudes)")
    ax.set_ylabel(r"wall-flux redistribution / value at $\phi=\pi/2$")
    ax.set_title(r"Phase selection of the quadrature pumping channel"
                 "\n" r"($k=0.8$, $Re=10$, $Pe=7$, $Gr=100$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_07_sin_phi.png")
    print("wrote FIG_v2_07_sin_phi.png in %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()

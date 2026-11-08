#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Campaign driver for the "physics-first" revision of
Symmetry-Selected Nonlinear Thermal Pumping in Periodically Heated
Mixed-Convection Channels.

Blocks implemented (perturbative part of the revision plan):
  D.1  Re_c(k) cartography        (7 k x 8 Re + Brent refinement)
  D.2  Gr-independence of Re_c    (4 Gr x Re scan + Brent per Gr)
  C.1/D.5  exact phase law        (quadrature coefficient Pi_perp, sin(phi) law)
  C.2  quantitative phase analysis (phase lags + term-wise decomposition)

Authoritative numerical source: validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py
(solve_bvp CPU reference).  No GPU approximation is used.
Every figure is produced from data computed by this script.
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
OUT.mkdir(exist_ok=True)
SOLVER = WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py"

_spec = importlib.util.spec_from_file_location("v4", SOLVER)
v4 = importlib.util.module_from_spec(_spec)
sys.modules["v4"] = v4
_spec.loader.exec_module(v4)

U0 = v4.U0
PE_REF = 7.0
GR_REF = 100.0
K_REF = 0.8


# ---------------------------------------------------------------------------
# Low level: first-order solves + the complex mixed kernel
# ---------------------------------------------------------------------------

def first_order(k, Re, Pe, Gr):
    """Return the four unit-amplitude first-order mode solves."""
    ts = v4.solve_thermal_bvp(k, Pe, True)
    ta = v4.solve_thermal_bvp(k, Pe, False)
    hs = v4.solve_hydro_bvp(k, Re, Gr, ts)
    ha = v4.solve_hydro_bvp(k, Re, Gr, ta)
    return ts, ta, hs, ha


def mixed_kernel(k, Re, Pe, Gr, n_eval=1001):
    """Complex mixed forcing kernel of the symmetric x antisymmetric interaction.

    With wall amplitudes eps_s, eps_a (complex) the physically averaged mixed
    forcing is E(y) = c*X(y) + d*Y(y), where c = Re[eps_a*conj(eps_s)],
    d = Im[eps_a*conj(eps_s)] and X, Y are the two real profiles returned here.
    c*X reproduces v4.cross_forcing when the amplitudes are real.
    """
    ts, ta, hs, ha = first_order(k, Re, Pe, Gr)
    y = np.linspace(-1.0, 1.0, n_eval)
    Ts, Ta = ts.theta(y), ta.theta(y)
    Ts_y, Ta_y = ts.theta_y(y), ta.theta_y(y)
    Us, Ua = hs.U(y), ha.U(y)
    Vs, Va = hs.V(y), ha.V(y)

    Z1 = Ua * np.conjugate(-1j * k * Ts) + Va * np.conjugate(Ts_y)
    Z2 = Us * np.conjugate(-1j * k * Ta) + Vs * np.conjugate(Ta_y)
    X = 0.5 * (np.real(Z1) + np.real(Z2))
    Yk = 0.5 * (np.imag(Z2) - np.imag(Z1))
    return y, X, Yk, dict(ts=ts, ta=ta, hs=hs, ha=ha, Z1=Z1, Z2=Z2)


def flux_functional(y, forcing, Pe):
    """Wall fluxes of the second-order mean thermal field driven by `forcing`."""
    mean = v4.solve_mean_bvp(y, forcing, Pe)
    thy = mean.theta_y(y)
    nu_plus = -thy[-1]
    nu_minus = thy[0]
    return 0.5 * (nu_plus - nu_minus), nu_plus, nu_minus, mean


def coefficients(k, Re, Pe, Gr, n_eval=1001):
    """In-phase Pi_T (published) and quadrature Pi_perp, plus consistency checks."""
    y, X, Yk, d = mixed_kernel(k, Re, Pe, Gr, n_eval)
    pi_c, nup_c, num_c, _ = flux_functional(y, X, Pe)
    pi_p, nup_p, num_p, _ = flux_functional(y, Yk, Pe)
    return dict(k=k, Re=Re, Pe=Pe, Pr=Pe / Re, Gr=Gr,
                Pi_c=pi_c, Pi_p=pi_p,
                Nu_plus_c=nup_c, Nu_minus_c=num_c,
                anti_sym_err_c=abs(nup_c + num_c),
                anti_sym_err_p=abs(nup_p + num_p),
                y=y, X=X, Y=Yk, modes=d)


def solve_pi(k, Re, Pe, Gr):
    """Published pumping coefficient at fixed Pe (Pr = Pe/Re)."""
    r = v4.solve_case(k, Re, Pe / Re, Gr)
    return float(r["Pi_flux"])


def locate_crossing(k, Pe, Gr, re_lo, re_hi):
    """Brent refinement of Re_c where Pi_T(k, Re) changes sign."""
    f = lambda R: solve_pi(k, R, Pe, Gr)
    return brentq(f, re_lo, re_hi, xtol=1e-8, rtol=1e-10, maxiter=200)


def save(df_or_obj, name):
    p = OUT / name
    if isinstance(df_or_obj, pd.DataFrame):
        df_or_obj.to_csv(p, index=False)
    else:
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(df_or_obj, fh, indent=2, default=float)
    print("   saved", p.name, flush=True)


# ---------------------------------------------------------------------------
# D.1 : Re_c(k) cartography
# ---------------------------------------------------------------------------

def campaign_D1():
    ks = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
    res = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
    rows = []
    t0 = time.time()
    for k in ks:
        for R in res:
            pi = solve_pi(k, R, PE_REF, GR_REF)
            rows.append(dict(k=k, Pe=PE_REF, Gr=GR_REF, Re=R, Pr=PE_REF / R,
                             Pi_T=pi, A=pi / GR_REF))
            print("D.1 k=%.2f Re=%4.1f  Pi_T=% .6e  A=% .6e  [%.0fs]"
                  % (k, R, pi, pi / GR_REF, time.time() - t0), flush=True)
    df = pd.DataFrame(rows)
    save(df, "D1_Re_grid.csv")

    # sign changes on the grid, then Brent refinement
    cross = []
    for k in ks:
        sub = df[df.k == k].sort_values("Re")
        sig = np.sign(sub.Pi_T.values)
        for i in range(len(sub) - 1):
            if sig[i] == 0 or sig[i + 1] == 0 or sig[i] * sig[i + 1] < 0:
                lo = float(sub.Re.values[i])
                hi = float(sub.Re.values[i + 1])
                if sig[i] == 0:
                    rec_c = lo
                elif sig[i + 1] == 0:
                    rec_c = hi
                else:
                    rec_c = locate_crossing(k, PE_REF, GR_REF, lo, hi)
                cross.append(dict(k=k, Re_c=rec_c, bracket_lo=lo, bracket_hi=hi,
                                  Pi_at_lo=float(sub.Pi_T.values[i]),
                                  Pi_at_hi=float(sub.Pi_T.values[i + 1])))
                print("   Re_c(k=%.2f) = %.6f" % (k, rec_c), flush=True)
                break
    cdf = pd.DataFrame(cross)
    save(cdf, "D1_Re_c.csv")
    fig_re_c_k(cdf, ks)
    fig_pi_vs_re(df, ks)
    return df, cdf


def fig_pi_vs_re(df, ks):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for k in ks:
        sub = df[df.k == k].sort_values("Re")
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
    print("   saved FIG_v2_01_piT_vs_Re_kartsch.png", flush=True)


def fig_re_c_k(cdf, ks):
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.plot(cdf.k, cdf.Re_c, "s-", color="tab:red")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$Re_c$  ($\Pi_T=0$, $Pe=7$, $Gr=100$)")
    for _, r in cdf.iterrows():
        ax.annotate("%.2f" % r.Re_c, (r.k, r.Re_c), textcoords="offset points",
                    xytext=(6, 4), fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_02_Re_c_vs_k.png", dpi=200)
    plt.close(fig)
    print("   saved FIG_v2_02_Re_c_vs_k.png", flush=True)


# ---------------------------------------------------------------------------
# D.2 : Gr-independence of the reversal point (Corollary 1)
# ---------------------------------------------------------------------------

def campaign_D2():
    grs = [25.0, 50.0, 100.0, 200.0]
    res = [25.0, 30.0, 35.0]
    rows = []
    for Gr in grs:
        for R in res:
            pi = solve_pi(K_REF, R, PE_REF, Gr)
            rows.append(dict(Gr=Gr, Re=R, Pi_T=pi, A=pi / Gr))
            print("D.2 Gr=%5.1f Re=%4.1f  Pi_T=% .6e  A=% .6e" % (Gr, R, pi, pi / Gr), flush=True)
    df = pd.DataFrame(rows)
    save(df, "D2_Gr_invariance.csv")

    # invariance of A across Gr at each Re
    inv = []
    for R in res:
        sub = df[df.Re == R]
        vals = sub.A.values
        inv.append(dict(Re=R, A_mean=float(vals.mean()),
                        max_rel_spread=float(np.max(np.abs(vals - vals.mean())) / abs(vals.mean()))))
    save(inv, "D2_A_invariance.json")

    # Re_c per Gr (Brent), expect identical roots
    cross = []
    for Gr in grs:
        lo, hi = 25.0, 35.0
        f = lambda R: solve_pi(K_REF, R, PE_REF, Gr)
        rc = locate_crossing(K_REF, PE_REF, Gr, lo, hi)
        cross.append(dict(Gr=Gr, Re_c=rc))
        print("   Re_c(Gr=%5.1f) = %.9f" % (Gr, rc), flush=True)
    cdf = pd.DataFrame(cross)
    save(cdf, "D2_Re_c_vs_Gr.csv")
    spread = float(np.ptp(cdf.Re_c.values))
    save(dict(Re_c_values={float(r.Gr): float(r.Re_c) for r in cdf.itertuples()},
              abs_spread=spread,
              rel_spread=spread / float(cdf.Re_c.mean())),
         "D2_Re_c_spread.json")

    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    for Gr in grs:
        sub = df[df.Gr == Gr]
        ax.plot(sub.Re, sub.A, "o-", label=r"$\mathrm{Gr}=%d$" % Gr)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"$Re$")
    ax.set_ylabel(r"$A=\Pi_T/\mathrm{Gr}$   ($k=0.8$, $Pe=7$)")
    ax.legend(fontsize=8)
    ax.set_title("Reversal point independent of Gr (Corollary 1)")
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_03_Gr_invariance.png", dpi=200)
    plt.close(fig)
    print("   saved FIG_v2_03_Gr_invariance.png", flush=True)
    return df, cdf


# ---------------------------------------------------------------------------
# C.1 / D.5 : exact phase law and quadrature coefficient
# ---------------------------------------------------------------------------

def analytic_coeffs(A_s, A_a):
    c = float(np.real(A_a * np.conjugate(A_s)))
    d = float(np.imag(A_a * np.conjugate(A_s)))
    return c, d


def campaign_phase():
    # reference case
    co = coefficients(K_REF, 10.0, PE_REF, GR_REF)
    ref = dict(k=K_REF, Re=10.0, Pe=PE_REF, Gr=GR_REF,
               Pi_c=co["Pi_c"], Pi_p=co["Pi_p"],
               check_against_solve_case=None,
               anti_sym_err_c=co["anti_sym_err_c"],
               anti_sym_err_p=co["anti_sym_err_p"])
    ref["check_against_solve_case"] = float(v4.solve_case(K_REF, 10.0, 0.7, GR_REF)["Pi_flux"])
    ref["relative_discrepancy"] = abs(ref["Pi_c"] - ref["check_against_solve_case"]) / abs(ref["check_against_solve_case"])
    print("PHASE ref: Pi_c=% .8e  Pi_perp=% .8e  ratio=% .6f  check=%.1e"
          % (ref["Pi_c"], ref["Pi_p"], ref["Pi_p"] / ref["Pi_c"], ref["relative_discrepancy"]), flush=True)

    # sin(phi) law, equal wall amplitudes eps=1
    eps = 1.0
    phis = [math.pi / 6, math.pi / 4, math.pi / 3, math.pi / 2,
            2 * math.pi / 3, 3 * math.pi / 4, 5 * math.pi / 6]
    rows = []
    for phi in phis:
        A_s = eps * math.cos(phi / 2) * np.exp(-1j * phi / 2)
        A_a = 1j * eps * math.sin(phi / 2) * np.exp(-1j * phi / 2)
        c, d = analytic_coeffs(A_s, A_a)
        pi_pred = c * ref["Pi_c"] + d * ref["Pi_p"]
        # direct assembly through the forcing profile (independent of the algebra)
        y, X, Yk = co["y"], co["X"], co["Y"]
        pi_direct, nup, num, _ = flux_functional(y, c * X + d * Yk, PE_REF)
        rows.append(dict(phi=phi, phi_over_pi=phi / math.pi, eps=eps,
                         Re_Aa_conj_As=c, Im_Aa_conj_As=d,
                         Pi_pred=pred_safe(pi_pred), Pi_direct=pi_direct,
                         abs_diff=abs(pi_pred - pi_direct),
                         sin_phi=math.sin(phi),
                         Pi_over_max=(eps**2 / 2) * pi_pred / ((eps**2 / 2) * ref["Pi_p"])))
        print("   phi=%.4f pi: c=%.3e d=%.3e  Pi=% .6e (direct % .6e)"
              % (phi / math.pi, c, d, pi_pred, pi_direct), flush=True)
    df = pd.DataFrame(rows)
    save(df, "D5_phase_law.csv")
    save(ref, "D5_reference_coefficients.json")

    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.plot(df.sin_phi, df.Pi_direct / df.Pi_direct.iloc[3], "o", label="computed")
    ax.plot([-1, 1], [-1, 1], "k-", lw=1, label=r"$y=x$ (exact identity)")
    ax.set_xlabel(r"$\sin\phi$")
    ax.set_ylabel(r"$\Pi_T(\phi)/\Pi_T(\pi/2)$")
    ax.legend(fontsize=8)
    ax.set_title(r"Exact phase selection, $k=0.8$, $Pe=7$, $Re=10$")
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_04_phase_law.png", dpi=200)
    plt.close(fig)

    # map Pi_perp / Pi_c over parameter space: when do they coincide?
    grid = []
    for k in [0.4, 0.8, 1.6, 3.2]:
        for R in [5.0, 10.0, 20.0, 40.0]:
            cc = coefficients(k, R, PE_REF, GR_REF)
            grid.append(dict(k=k, Re=R, Pe=PE_REF, Gr=GR_REF,
                             Pi_c=cc["Pi_c"], Pi_p=cc["Pi_p"],
                             ratio=cc["Pi_p"] / cc["Pi_c"]))
            print("   ratio k=%.1f Re=%4.1f : %.6f" % (k, R, grid[-1]["ratio"]), flush=True)
    save(pd.DataFrame(grid), "D5_ratio_map.csv")
    return df, ref


def pred_safe(x):
    return float(x)


# ---------------------------------------------------------------------------
# C.2 : quantitative phase analysis of the Re-driven reversal
# ---------------------------------------------------------------------------

def campaign_phase_lag():
    res = [5.0, 8.0, 10.0, 15.0, 20.0, 25.0, 28.0, 30.0, 32.0, 35.0, 40.0]
    rows = []
    for R in res:
        y, X, Yk, d = mixed_kernel(K_REF, R, PE_REF, GR_REF)
        ts, ta, hs, ha = d["ts"], d["ta"], d["hs"], d["ha"]
        yc = np.array([0.0, 0.5])
        th_s, th_a = ts.theta(yc), ta.theta(yc)
        v_s, v_a = hs.V(yc), ha.V(yc)
        pi_c, _, _, _ = flux_functional(y, X, PE_REF)
        # term-wise contributions to the published coefficient
        yv = y
        Ts, Ta = ts.theta(yv), ta.theta(yv)
        Ts_y, Ta_y = ts.theta_y(yv), ta.theta_y(yv)
        Us, Ua = hs.U(yv), ha.U(yv)
        Vs, Va = hs.V(yv), ha.V(yv)
        terms = {
            "Ua_dX_Ths": 0.5 * np.real(Ua * np.conjugate(-1j * K_REF * Ts)),
            "Us_dX_Tha": 0.5 * np.real(Us * np.conjugate(-1j * K_REF * Ta)),
            "Va_dY_Ths": 0.5 * np.real(Va * np.conjugate(Ts_y)),
            "Vs_dY_Tha": 0.5 * np.real(Vs * np.conjugate(Ta_y)),
        }
        contrib = {}
        for name, prof in terms.items():
            val, _, _, _ = flux_functional(yv, prof, PE_REF)
            contrib[name] = val
        row = dict(Re=R, k=K_REF, Pe=PE_REF, Gr=GR_REF, Pi_T=pi_c,
                   abs_Th_s_0=abs(th_s[0]), arg_Th_s_0=float(np.angle(th_s[0])),
                   abs_Th_a_0=abs(th_a[0]), arg_Th_a_0=float(np.angle(th_a[0])),
                   abs_V_s_0=abs(v_s[0]), arg_V_s_0=float(np.angle(v_s[0])),
                   abs_V_a_0=abs(v_a[0]), arg_V_a_0=float(np.angle(v_a[0])),
                   dphi_s=float(np.angle(v_s[0]) - np.angle(th_s[0])),
                   dphi_a=float(np.angle(v_a[0]) - np.angle(th_a[0])),
                   dphi_s_half=float(np.angle(v_s[1]) - np.angle(th_s[1])),
                   dphi_a_half=float(np.angle(v_a[1]) - np.angle(th_a[1])),
                   **contrib)
        rows.append(row)
        print("D.C2 Re=%4.1f  Pi_T=% .6e  dphi_s=%+ .4f dphi_a=%+ .4f  terms=%s"
              % (R, pi_c, row["dphi_s"], row["dphi_a"],
                 " ".join("%s=%.3e" % (kk, vv) for kk, vv in contrib.items())), flush=True)
    df = pd.DataFrame(rows)
    save(df, "C2_phase_lag.csv")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax = axes[0]
    ax.plot(df.Re, df.Pi_T, "o-", color="tab:blue", label=r"$\Pi_T$")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_xlabel(r"$Re$ ($k=0.8$, $Pe=7$)")
    ax.set_ylabel(r"$\Pi_T$", color="tab:blue")
    ax2 = ax.twinx()
    ax2.plot(df.Re, df.dphi_s, "s-", color="tab:red", label=r"$\Delta\varphi^{(s)}$")
    ax2.plot(df.Re, df.dphi_a, "^-", color="tab:green", label=r"$\Delta\varphi^{(a)}$")
    ax2.axhline(0.0, color="grey", ls="--", lw=0.8)
    ax2.set_ylabel(r"phase lag $\Delta\varphi$ [rad]", color="k")
    ax.set_title("Reversal vs velocity-temperature phase lag")
    ax = axes[1]
    for name, col in [("Ua_dX_Ths", "tab:red"), ("Us_dX_Tha", "tab:orange"),
                      ("Va_dY_Ths", "tab:purple"), ("Vs_dY_Tha", "tab:brown")]:
        ax.plot(df.Re, df[name], "o-", ms=4, color=col, label=name.replace("_", " "))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"$Re$")
    ax.set_ylabel(r"contribution to $\Pi_T$")
    ax.legend(fontsize=7)
    ax.set_title("Term-wise decomposition of the mixed forcing")
    fig.tight_layout()
    fig.savefig(OUT / "FIG_v2_05_phase_analysis.png", dpi=200)
    plt.close(fig)
    print("   saved FIG_v2_05_phase_analysis.png", flush=True)
    return df


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if which in ("all", "D1"):
        campaign_D1()
    if which in ("all", "D2"):
        campaign_D2()
    if which in ("all", "phase"):
        campaign_phase()
    if which in ("all", "C2"):
        campaign_phase_lag()
    print("TOTAL %.0f s" % (time.time() - t0), flush=True)

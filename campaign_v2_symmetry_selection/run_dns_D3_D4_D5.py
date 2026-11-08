#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepared-but-not-launched DNS campaigns D.3, D.4 and D.5 of the action plan.

The case matrix is not the literal grid of the plan (k = 0.6, 0.8, 1.2 x
Re = 20, 30, 35): that grid provably fails its own purpose, because the
perturbative cartography gives Re_c(0.6) = 51.25 and Re_c(1.2) = 15.09, so all
nine of its cells sit on a single side of the crossing at k = 0.6 and k = 1.2 and
would measure nothing about a reversal. It is also unresolvable near the
crossing. Re-analysing the published fixed-Pe DNS sweep against the exact
perturbative coefficients shows a deviation that is flat at about 2.7e-4 in
Pi_T units for small coefficients (which is why the stored Re = 30 point returns
the wrong sign: +7.5e-7 measured against -2.0e-4 predicted) and proportional to
about 1.5 % of Pi_T for large ones. Every case therefore carries two resolvability
margins: H1 assumes the deviation is fixed in coefficient space, H2 assumes it is
absolute in q_+ - q_- and hence improves with the forcing amplitude. Cells are
kept when the sign is resolvable under H1; the single low-amplitude ladder run
(eps_a = 0.02 at Re = 30, where H1 says no and H2 says yes) is deliberately kept
because it decides between the two models in either direction.

    D.3  reversal bracketed at three wavenumbers, both sides of Re_c   -> 6 runs
    D.3b low-k bracket at k = 0.6 (positive branch capped at 6.6e-4)   -> 2 runs
    D.3c near-crossing amplitude ladder, error-model discriminator     -> 1 run
    D.4  control point far from the crossing (Re = 50, k = 0.8)        -> 1 run
    D.5  quadrature test: eps_1 = eps_2 = 0.03, phi = 0, +-pi/2        -> 3 runs

Pe = 7 and Gr = 100 throughout, Pr = Pe/Re (fixed-Pe convention of the
manuscript), mesh and time step of the validated reference runs.

Every run is checkpointed as soon as it finishes and skipped on a later
invocation, so an interrupted campaign resumes instead of restarting. The script
is a DRY RUN unless --go is given; the dry run only solves the perturbative BVP
chain (seconds) and prints the matrix and the falsifiable predictions.

Launch (only after explicit approval):
    py -X utf8 campaign_v2_symmetry_selection/run_dns_D3_D4_D5.py --go
Single case:
    py -X utf8 campaign_v2_symmetry_selection/run_dns_D3_D4_D5.py --go --only D3_k1.0_Re15
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path(r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification")
OUT = WORK / "campaign_v2_symmetry_selection"
CHECKPOINT = OUT / "DNS_D345_results.csv"
PREDICT = OUT / "DNS_D345_predictions.csv"
LOWK_CSV = OUT / "D1_Re_c_lowk_extension.csv"

sys.path.insert(0, str(WORK))
_spec = importlib.util.spec_from_file_location(
    "v4", WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py")
v4 = importlib.util.module_from_spec(_spec)
sys.modules["v4"] = v4
_spec.loader.exec_module(v4)

PE, GR = 7.0, 100.0
EPS_S = 0.03
EPS_A = 0.005
AMP_REF = EPS_A * EPS_S  # bilinear amplitude product of the validated reference runs
# DNS-vs-perturbative deviation calibrated on dns_re_sweep.csv: the deviation is
# flat at about 2.7e-4 in Pi_T units for small coefficients and grows by about
# 1.5 % of |Pi_T| for large ones.
DNS_DEV = 2.71e-4
DNS_DEV_REL = 1.52e-2
# mesh and time step of the already-validated reference runs
NX, NY, DT, STEPS = 64, 49, 0.0015, 30000


def dns_deviation(pi_t):
    """Predicted absolute error on a recovered pump coefficient."""
    return max(DNS_DEV, DNS_DEV_REL * abs(pi_t))


# -----------------------------------------------------------------------------
# perturbative reference (DNS is imported only under --go: the dry run stays cheap)
# -----------------------------------------------------------------------------
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


_COEFF_CACHE = {}


def coefficients(k, Re, Pe=PE, Gr=GR):
    """(Pi_c, Pi_p) and the antisymmetry residuals at one parameter point."""
    key = (round(k, 6), round(Re, 6), round(Pe, 6), round(Gr, 6))
    if key not in _COEFF_CACHE:
        y, Kc, Kq = kernels(k, Re, Pe, Gr)
        pi_c, anti_c = wall_flux(y, Kc, Pe)
        pi_p, anti_p = wall_flux(y, Kq, Pe)
        _COEFF_CACHE[key] = (pi_c, pi_p, anti_c, anti_p)
    return _COEFF_CACHE[key]


def calibrate_floor():
    """|Pi_T^DNS - Pi_T^pert| over the published fixed-Pe Re sweep at k = 0.8."""
    path = WORK / "dns_physical_validation_suite" / "dns_re_sweep.csv"
    if not path.exists():
        return None
    rows = []
    for _, r in pd.read_csv(path).iterrows():
        pi_c = coefficients(0.8, float(r.Re))[0]
        rows.append((float(r.Re), pi_c, float(r.PiT) - pi_c))
    return rows


def bilinear_factors(eps1, eps2, phi):
    """c = Re(A_a A_s*), d = Im(A_a A_s*) with theta_hat(-1)=eps2 e^{-i phi}."""
    A_s = 0.5 * (eps1 + eps2 * np.exp(-1j * phi))
    A_a = 0.5 * (eps1 - eps2 * np.exp(-1j * phi))
    p = A_a * np.conjugate(A_s)
    return float(np.real(p)), float(np.imag(p))


def expected_row(case):
    c, d = bilinear_factors(case["eps1"], case["eps2"], case["phi"])
    pi_c, pi_p, anti_c, anti_p = coefficients(case["k"], case["Re"])
    # two-coefficient phase law: Nu_2,0(+1) = c Pi_c + d Pi_p, q_+ - q_- = 2 Nu_2,0(+1)
    q_pred = 2.0 * (c * pi_c + d * pi_p)
    amp = math.hypot(c, d)
    proj = q_pred / (2.0 * amp) if amp > 0 else np.nan
    # H1: the calibrated deviation lives in coefficient space (amplitude-independent)
    # H2: it is absolute in q_+ - q_- and was measured at the reference amplitude
    m_h1 = abs(proj) / dns_deviation(proj) if np.isfinite(proj) else np.nan
    m_h2 = m_h1 * (amp / AMP_REF) if np.isfinite(m_h1) else np.nan
    return dict(c=c, d=d, PiT_perturbative=pi_c, PiT_perp=pi_p,
                q_diff_predicted=q_pred, PiT_projected=proj,
                margin_H1=m_h1, margin_H2=m_h2,
                expected_sign=("zero" if abs(q_pred) < 1e-14 else
                               ("positive" if q_pred > 0 else "negative")),
                antisym_err=max(anti_c, anti_p))


# -----------------------------------------------------------------------------
# case matrix
# -----------------------------------------------------------------------------
def mk(case_id, family, k, Re, phi=0.0, eps_a=EPS_A, purpose=""):
    return dict(case_id=case_id, family=family, purpose=purpose, k=k, Re=Re,
                Pr=PE / Re, phi=phi, eps_a=eps_a,
                eps1=EPS_S + eps_a, eps2=EPS_S - eps_a,
                Gr=GR, Pe=PE, Nx=NX, Ny=NY, dt=DT, steps=STEPS)


def build_cases():
    cases = [
        # D.3: bracket of the crossing at three wavenumbers (Re_c from D.1)
        mk("D3_k0.8_Re20", "D.3", 0.8, 20.0, purpose="below Re_c=30.838"),
        mk("D3_k0.8_Re40", "D.3", 0.8, 40.0, purpose="above Re_c=30.838"),
        mk("D3_k1.0_Re15", "D.3", 1.0, 15.0, purpose="below Re_c=20.487"),
        mk("D3_k1.0_Re25", "D.3", 1.0, 25.0, purpose="above Re_c=20.487"),
        mk("D3_k1.2_Re10", "D.3", 1.2, 10.0, purpose="below Re_c=15.089"),
        mk("D3_k1.2_Re20", "D.3", 1.2, 20.0, purpose="above Re_c=15.089"),
        # D.3b: low-k bracket, positive branch capped at 6.6e-4 -> marginal
        mk("D3b_k0.6_Re40", "D.3b", 0.6, 40.0, purpose="below Re_c=51.252"),
        mk("D3b_k0.6_Re70", "D.3b", 0.6, 70.0, purpose="above Re_c=51.252"),
        # D.3c: amplitude ladder at the crossing, discriminates the error model
        mk("D3c_k0.8_Re30_ea02", "D.3c", 0.8, 30.0, eps_a=0.02,
           purpose="near-crossing ladder: H1 says unresolvable, H2 says yes"),
        # D.4: control far from the crossing
        mk("D4_k0.8_Re50", "D.4", 0.8, 50.0, purpose="far above Re_c(0.8)"),
    ]
    for phi, tag in ((0.0, "0"), (math.pi / 2, "plus_pi_over_2"),
                     (-math.pi / 2, "minus_pi_over_2")):
        # equal wall amplitudes with a phase lag: pure quadrature forcing,
        # eps_a = (eps1-eps2)/2 = 0, so the run_dns estimator is undefined by
        # construction and only q_+ - q_- is used (normalised below)
        cases.append(dict(case_id="D5_phi_%s" % tag, family="D.5",
                          purpose="sin(phi) selection rule", k=0.8, Re=10.0,
                          Pr=0.7, phi=phi, eps_a=0.0, eps1=0.03, eps2=0.03,
                          Gr=GR, Pe=PE, Nx=NX, Ny=NY, dt=DT, steps=STEPS))
    return cases


def refine_low_k_re_c():
    """Re_c at the two wavenumbers the D.1 Re grid (<= 40) could not bracket."""
    from scipy.optimize import brentq
    rows = []
    for k, lo, hi in ((0.6, 45.0, 90.0), (0.4, 60.0, 120.0)):
        f = lambda R, k=k: coefficients(k, R)[0]
        rc = brentq(f, lo, hi, xtol=1e-8, rtol=1e-10, maxiter=200)
        rows.append(dict(k=k, Re_c=rc, bracket_lo=lo, bracket_hi=hi,
                         Pe=PE, Gr=GR, Pr_at_Re_c=PE / rc))
        print("   Re_c(k=%.2f) = %.6f  (Pr = %.4f)" % (k, rc, PE / rc))
    pd.DataFrame(rows).to_csv(LOWK_CSV, index=False)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--go", action="store_true",
                    help="actually integrate (default is a dry run)")
    ap.add_argument("--only", default=None, help="run a single case_id")
    ap.add_argument("--force", action="store_true",
                    help="redo cases already present in the checkpoint")
    ap.add_argument("--lowk", action="store_true",
                    help="also write the low-k Re_c extension of the D.1 map")
    args = ap.parse_args()

    t_start = time.time()
    cal = calibrate_floor()
    if cal:
        print("error model, published DNS Re sweep at k = 0.8 (eps_a = 0.005):")
        for Re, pi_c, e in cal:
            print("   Re=%-5.1f |Pi_T^pert|=%+.4e  Pi_T^DNS - Pi_T^pert = %+.4e"
                  % (Re, pi_c, e))
        flat = [abs(e) for _, p, e in cal if abs(p) < 1e-3]
        rel = [abs(e / p) for _, p, e in cal if abs(p) > 1e-2]
        print("   derived flat deviation %.3e (constants %.3e), relative %.3e (%.3e)"
              % (max(flat), DNS_DEV, max(rel), DNS_DEV_REL))
    print("   H1 = deviation fixed in coefficient space, "
          "H2 = deviation absolute in q_+-q_- at the reference amplitude")

    cases = build_cases()
    if args.only:
        cases = [c for c in cases if c["case_id"] == args.only]
        if not cases:
            sys.exit("unknown case_id %r" % args.only)

    done = set()
    if CHECKPOINT.exists() and not args.force:
        done = set(pd.read_csv(CHECKPOINT).case_id)

    pi_c0, pi_p0, _, _ = coefficients(0.8, 10.0)
    ref = json.load(open(OUT / "D5_reference_coefficients.json"))
    print("\nreference coefficients (k=0.8, Pe=7, Re=10, Pr=0.7, Gr=100):")
    print("  Pi_T       = %+.8e   (stored %+.8e, rel diff %.2e)"
          % (pi_c0, ref["Pi_c"], abs(pi_c0 - ref["Pi_c"]) / abs(ref["Pi_c"])))
    print("  Pi_T_perp  = %+.8e   (stored %+.8e, rel diff %.2e)"
          % (pi_p0, ref["Pi_p"], abs(pi_p0 - ref["Pi_p"]) / abs(ref["Pi_p"])))

    print("\nD.5 falsifiable predictions (q_+ - q_- = 2 d Pi_T_perp):")
    for phi in (0.0, math.pi / 2, -math.pi / 2):
        _, d = bilinear_factors(0.03, 0.03, phi)
        print("  phi=%+6.3f  d=%.4e  q_+-q_- predicted %+.6e"
              % (phi, d, 2.0 * d * pi_p0))

    if args.lowk:
        print("\nlow-k refinement of the Re_c cartography:")
        refine_low_k_re_c()

    print("\ncase matrix (%d cases, %d already checkpointed)"
          % (len(cases), len([c for c in cases if c["case_id"] in done])))
    hdr = ("%-20s %-4s %-4s %-6s %-7s %-7s %-8s %-12s %-9s %-7s %-7s %s"
           % ("case_id", "fam", "k", "Re", "Pr", "phi", "eps_a", "q_pred",
              "expected", "margH1", "margH2", "asym.err"))
    print(hdr)
    print("-" * len(hdr))
    pending = []
    all_exp = []
    for c in cases:
        exp = expected_row(c)
        all_exp.append(exp)
        print("%-20s %-4s %-4.1f %-6.1f %-7.4f %-7.4f %-8.4f %-12s %-9s "
              "%-7.1f %-7.1f %.1e"
              % (c["case_id"], c["family"], c["k"], c["Re"], c["Pr"], c["phi"],
                 c["eps_a"], "%.4e" % exp["q_diff_predicted"],
                 exp["expected_sign"], exp["margin_H1"], exp["margin_H2"],
                 exp["antisym_err"]))
        if c["case_id"] not in done:
            pending.append((c, exp))

    dead = [c["case_id"] for c, e in pending
            if np.isfinite(e["margin_H2"]) and max(e["margin_H1"],
                                                   e["margin_H2"]) < 2.0]
    disc = [c["case_id"] for c, e in pending
            if np.isfinite(e["margin_H1"]) and e["margin_H1"] < 2.0
            and e["margin_H2"] >= 2.0]
    if dead:
        print("NOTE unresolvable under both error models: %s" % ", ".join(dead))
    if disc:
        print("NOTE error-model discriminators (H1 no, H2 yes): %s"
              % ", ".join(disc))

    pred_rows = []
    for c, exp in zip(cases, all_exp):
        pred_rows.append(dict(case_id=c["case_id"], family=c["family"],
                              purpose=c["purpose"], k=c["k"], Re=c["Re"],
                              Pr=c["Pr"], Pe=c["Pe"], Gr=c["Gr"],
                              phi=c["phi"], eps1=c["eps1"], eps2=c["eps2"],
                              eps_a=c["eps_a"], c_bilinear=exp["c"],
                              d_bilinear=exp["d"],
                              PiT_perturbative=exp["PiT_perturbative"],
                              PiT_perp_perturbative=exp["PiT_perp"],
                              PiT_projected=exp["PiT_projected"],
                              q_diff_predicted=exp["q_diff_predicted"],
                              expected_sign=exp["expected_sign"],
                              margin_H1=exp["margin_H1"],
                              margin_H2=exp["margin_H2"],
                              antisym_err=exp["antisym_err"],
                              checkpointed=c["case_id"] in done))
    pd.DataFrame(pred_rows).to_csv(PREDICT, index=False)
    print("\npredictions written to", PREDICT.name)

    print("\nperturbative references computed in %.1f s" % (time.time() - t_start))
    print("estimated cost: %d runs x %d steps at Nx=%d, Ny=%d, dt=%g"
          % (len(pending), STEPS, NX, NY, DT))
    for c, exp in pending:
        print("  %-20s %s | predicted q_+-q_- %+.6e (%s)"
              % (c["case_id"], c["purpose"], exp["q_diff_predicted"],
                 exp["expected_sign"]))
    if not args.go:
        print("DRY RUN: nothing integrated. Re-invoke with --go to launch.")
        return 0

    from dns_physical_validation_suite import run_dns  # noqa: E402

    new_file = not CHECKPOINT.exists()
    handle = open(CHECKPOINT, "a", newline="")
    for c, exp in pending:
        t0 = time.time()
        res = run_dns(eps1=c["eps1"], eps2=c["eps2"], k=c["k"], Re=c["Re"],
                      Pr=c["Pr"], Gr=c["Gr"], Nx=c["Nx"], Ny=c["Ny"],
                      dt=c["dt"], steps=c["steps"], phi=c["phi"],
                      verbose=False)
        q_diff = res["q_top"] - res["q_bottom"]
        amp = math.hypot(exp["c"], exp["d"])
        row = dict(case_id=c["case_id"], family=c["family"],
                   purpose=c["purpose"], k=c["k"], Re=c["Re"], Pr=c["Pr"],
                   Pe=c["Pe"], Gr=c["Gr"], phi=c["phi"], eps1=c["eps1"],
                   eps2=c["eps2"], eps_a=c["eps_a"], Nx=c["Nx"], Ny=c["Ny"],
                   dt=c["dt"], steps=c["steps"],
                   c_bilinear=exp["c"], d_bilinear=exp["d"],
                   PiT_perturbative=exp["PiT_perturbative"],
                   PiT_perp_perturbative=exp["PiT_perp"],
                   PiT_projected_predicted=exp["PiT_projected"],
                   q_diff_predicted=exp["q_diff_predicted"],
                   expected_sign=exp["expected_sign"],
                   margin_H1=exp["margin_H1"], margin_H2=exp["margin_H2"],
                   q_top=res["q_top"], q_bottom=res["q_bottom"],
                   q_diff_dns=q_diff,
                   PiT_dns_projected=(q_diff / (2.0 * amp) if amp > 0 else np.nan),
                   PiT_dns_suite=res["PiT"],
                   final_residual=res["final_residual"],
                   wall_seconds=time.time() - t0)
        pd.DataFrame([row]).to_csv(handle, mode="a",
                                   header=new_file, index=False)
        new_file = False
        handle.flush()
        os.fsync(handle.fileno())
        print("DONE %-20s q_diff=%+.6e (predicted %+.6e) resid=%.2e %.0fs"
              % (c["case_id"], q_diff, exp["q_diff_predicted"],
                 res["final_residual"], time.time() - t0), flush=True)
    handle.close()
    print("\ncheckpoint:", CHECKPOINT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

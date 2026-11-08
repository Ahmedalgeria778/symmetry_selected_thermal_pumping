"""Cross-check every number quoted in main_v2.tex against the stored data.

Run from the project root:
    py -X utf8 campaign_v2_symmetry_selection/verify_manuscript_numbers.py
Exits non-zero if any quoted value disagrees with the data files.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DENSE = os.path.join(ROOT, "validation_outputs_ULTIMATE_DENSE_CPU")

fails = []


def check(label, got, want, tol, rel=False):
    ok = abs(got - want) <= (tol * abs(want) if rel else tol)
    if not ok:
        fails.append(label)
    print(("PASS  " if ok else "FAIL  ") + f"{label}: {got!r} vs {want!r}")


def true(label, got, want=True):
    ok = bool(got is want or got == want)
    if not ok:
        fails.append(label)
    print(("PASS  " if ok else "FAIL  ") + f"{label}: {got!r}")


# ------------------------------------------------------ reference coefficients
ref = json.load(open(os.path.join(HERE, "D5_reference_coefficients.json")))
piT, piTp = ref["Pi_c"], ref["Pi_p"]
check("PiT = -4.74850846e-2", piT, -4.74850846e-2, 1e-8, rel=True)
check("PiT_perp = +8.89496095e-2", piTp, +8.89496095e-2, 1e-8, rel=True)
check("ratio = -1.8732", piTp / piT, -1.8732, 5e-5)
check("kernel vs published solver = 4.3e-13", ref["relative_discrepancy"],
      4.3e-13, 5e-14)
true("antisymmetry residual (in-phase) <= 6.9e-18",
     ref["anti_sym_err_c"] <= 6.94e-18)
true("antisymmetry residual (quadrature) <= 1.4e-17",
     ref["anti_sym_err_p"] <= 1.39e-17)

# ------------------------------------------------------------- DNS predictions
check("Nu(+1) equal-amplitude phi=pi/2 = +4.003e-5",
      0.5 * 0.03 ** 2 * piTp, 4.003e-5, 5e-9)
check("Nu(+1) eps_s*eps_a = -7.123e-6", 0.03 * 0.005 * piT, -7.123e-6, 5e-10)
check("magnitude factor 5.6",
      abs(0.5 * 0.03 ** 2 * piTp / (0.03 * 0.005 * piT)), 5.62, 0.02)

# ------------------------------------------------------------- D1 Re_c(k) map
d1 = pd.read_csv(os.path.join(HERE, "D1_Re_c.csv"))
grid = pd.read_csv(os.path.join(HERE, "D1_Re_grid.csv"))
true("D1 covers 7 wavenumbers x 8 Reynolds", grid.shape[0], 56)
for k, want in [(0.8, 30.838497), (1.0, 20.487029), (1.2, 15.089038),
                (1.5, 11.069430), (2.0, 8.541316)]:
    got = float(d1.loc[d1.k == k, "Re_c"].iloc[0])
    check(f"Re_c(k={k})", got, want, 5e-7)
for k in (0.4, 0.6):
    sub = grid[grid.k == k]
    true(f"no sign change at k={k} over Re in [5,40]",
         bool((sub.Pi_T < 0).all()))
kk, rr = d1.k.values, d1.Re_c.values
a, b = np.polyfit(np.log(kk), np.log(rr), 1)
dev = np.abs(np.exp(b) * kk ** a - rr) / rr * 100.0
check("fit prefactor 20.91", np.exp(b), 20.91, 0.01)
check("fit exponent -1.41", a, -1.41, 0.01)
check("fit max deviation 7.7 %", dev.max(), 7.7, 0.1)

# --------------------------------------------------------- D2 Gr invariance
sp = json.load(open(os.path.join(HERE, "D2_Re_c_spread.json")))
ainv = json.load(open(os.path.join(HERE, "D2_A_invariance.json")))
check("Re_c absolute spread 1.2e-10", sp["abs_spread"], 1.234e-10, 5e-13)
check("Re_c value 30.838496520", max(sp["Re_c_values"].values()),
      30.838496520, 5e-9)
check("worst-case A relative spread 7.3e-11",
      max(e["max_rel_spread"] for e in ainv), 7.291e-11, 5e-14)
check("A(Re=25) = -2.089882e-5",
      [e for e in ainv if e["Re"] == 25][0]["A_mean"], -2.089882e-5, 1e-11)
check("A(Re=30) = -1.958845e-6",
      [e for e in ainv if e["Re"] == 30][0]["A_mean"], -1.958845e-6, 1e-12)
check("A(Re=35) = +6.745165e-6",
      [e for e in ainv if e["Re"] == 35][0]["A_mean"], +6.745165e-6, 1e-12)
gr = pd.read_csv(os.path.join(HERE, "P2_Gr_linearity_two_coeffs.csv"))
for col, want in [("A_c", 1.02e-13), ("A_p", 3.68e-13)]:
    v = gr[col].values
    check(f"{col} relative spread over Gr=25-400",
          (v.max() - v.min()) / abs(v.mean()), want, 5e-15)

# ----------------------------------------------- two-branch Re scan and ratio
scan = pd.read_csv(os.path.join(HERE, "P2_Re_scan_two_coeffs.csv"))
check("PiT_perp(Re=5) = 1.86e-1", scan.Pi_p.iloc[0], 1.86e-1, 2e-3)
check("PiT_perp(Re=50) = 8.40e-3", scan.Pi_p.iloc[-1], 8.40e-3, 2e-4)
true("PiT_perp > 0 for all scanned Re", bool((scan.Pi_p > 0).all()))
true("PiT_perp decreasing over the scan", bool(np.all(np.diff(scan.Pi_p) < 0)))
true("in-phase branch crosses once", bool((np.diff(np.sign(scan.Pi_c)) != 0).sum() == 1))
branches = json.load(open(os.path.join(HERE, "P2_Re_c_two_branches.json")))
check("in-phase Re_c = 30.8385", branches["Re_c_inphase"], 30.838497, 5e-7)
true("quadrature branch has no Re_c", branches["Re_c_quadrature"] is None)
ratio = pd.read_csv(os.path.join(HERE, "D5_ratio_map.csv"))
check("ratio min -36.3", ratio.ratio.min(), -36.264, 0.01)
check("ratio max +11.0", ratio.ratio.max(), 11.023, 0.01)
check("ratio at (0.4,5) = -1.17",
      ratio[(ratio.k == 0.4) & (ratio.Re == 5)].ratio.iloc[0], -1.171, 0.002)

# ----------------------------------------------------- two-coefficient law test
val = pd.read_csv(os.path.join(HERE, "P2_two_coeff_law_validation.csv"))
true("two-coefficient law tested on 12 randomised triples", val.shape[0], 12)
check("two-coefficient law max relative deviation 5.0e-12",
      val.rel_err.abs().max(), 4.971e-12, 1e-14)
check("two-coefficient law max absolute deviation 1.5e-14",
      val.abs_err.abs().max(), 1.473e-14, 1e-16)
law = pd.read_csv(os.path.join(HERE, "D5_phase_sweep_direct.csv"))
check("phase sweep covers phi in [0, pi] (lower)", law.phi_over_pi.min(), 0.0, 1e-12)
check("phase sweep covers phi in [0, pi] (upper)", law.phi_over_pi.max(), 1.0, 1e-12)
check("phase sweep has 25 points", law.shape[0], 25, 0)
check("sin(phi) law max deviation 2.4e-13",
      law.dev_from_sin.abs().max(), 2.367e-13, 1e-15)
check("direct vs algebraic max abs error 1.1e-14",
      law.abs_err_pred.max(), 1.100e-14, 1e-16)
check("in-phase invariant Re(A_a As*) vanishes to 5.6e-17",
      law.Re_Aa_As.abs().max(), 5.551e-17, 1e-18)
true("sweep antisymmetry residual <= 1.4e-17", law.antisym_err.max() <= 1.39e-17)
check("quadrature invariant = (eps^2/2) sin(phi)",
      (law.Im_Aa_As - 0.5 * law.eps ** 2 * law.sin_phi).abs().max(), 0.0, 1e-15)

# ------------------------------------------------------------- phase analysis
unw = pd.read_csv(os.path.join(HERE, "C2_phase_unwrapped.csv"))
check("term-balance max absolute error 1.4e-14",
      unw.term_balance_err.abs().max(), 1.3684e-14, 1e-17)
check("term-balance max relative error 1.6e-5",
      (unw.term_balance_err / unw.Pi_T).abs().max(), 1.59e-5, 1e-6)
row5, row50 = unw.iloc[0], unw.iloc[-1]
check("|T_Va|/|T_Ua| at Re=5 = 3.40",
      abs(row5.T_Va_dYs / row5.T_Ua_dXs), 3.396, 0.01)
check("|T_Va|/|T_Ua| at Re=50 = 0.936",
      abs(row50.T_Va_dYs / row50.T_Ua_dXs), 0.9361, 0.002)
rat = (unw.T_Va_dYs.abs() / unw.T_Ua_dXs.abs()).values
true("dominant-term ratio decreases monotonically", bool(np.all(np.diff(rat) < 0)))
true("T_Va_dYs never changes sign", bool((unw.T_Va_dYs < 0).all()))
true("T_Ua_dXs never changes sign", bool((unw.T_Ua_dXs > 0).all()))
check("unwrapped dphi_a(0) start -2.154", unw.dphi_a_0_unw.iloc[0], -2.154, 5e-4)
check("unwrapped dphi_a(0) end -1.532", unw.dphi_a_0_unw.iloc[-1], -1.532, 5e-4)
true("unwrapped dphi_a(0) increases monotonically",
     bool(np.all(np.diff(unw.dphi_a_0_unw) > 0)))
check("sweep overlap start 0.993", unw.ov_a_dYs_0.iloc[0], 0.9925, 5e-4)
check("sweep overlap end 0.878", unw.ov_a_dYs_0.iloc[-1], 0.8781, 5e-4)
rc = unw.iloc[(unw.Re - 30.838497).abs().argsort().iloc[0]]
check("subdominant net at crossing +2.7e-4", rc.T_Us_dXa + rc.T_Vs_dYa,
      2.669e-4, 1e-6)
check("two-term proxy crossing Re=32.713",
      json.load(open(os.path.join(HERE, "C2_phase_crossings.json")))
      ["term_T_Va_dYs_plus_T_Ua_zero"][0], 32.713, 5e-4)
for re_, want in [(5, (-1.5951e-1, 5.504e-2, -1.870e-1, -2.537e-2, -2.236e-3)),
                  (10, (-4.7485e-2, 5.899e-2, -9.726e-2, -2.707e-3, -6.504e-3)),
                  (20, (-6.428e-3, 3.737e-2, -4.251e-2, 4.895e-3, -6.180e-3)),
                  (30, (-1.959e-4, 2.378e-2, -2.420e-2, 4.829e-3, -4.611e-3)),
                  (40, (1.081e-3, 1.655e-2, -1.599e-2, 3.939e-3, -3.420e-3)),
                  (50, (1.330e-3, 1.232e-2, -1.154e-2, 3.157e-3, -2.614e-3))]:
    r = unw.iloc[(unw.Re - re_).abs().argsort().iloc[0]]
    got = (r.Pi_T, r.T_Ua_dXs, r.T_Va_dYs, r.T_Us_dXa, r.T_Vs_dYa)
    for name, g, w in zip(["PiT", "T_Ua", "T_Va", "T_Us", "T_Vs"], got, want):
        check(f"tab:term_balance Re={re_} {name}", g, w, 0.5e-3 * abs(w) + 1e-9)

# -------------------------------------------------------------- Pe_c endpoint
dense = pd.read_csv(os.path.join(DENSE, "ultimate_dense_cpu_scan.csv"))
sub = dense[np.isclose(dense["Pe"], 10.9)]
check("PiT(10,10.9) = -1.4068443e-6",
      float(sub.loc[sub.k.idxmax(), "Pi_T"]), -1.4068443e-6, 1e-7, rel=True)
crit = dense[np.isclose(dense["Pe"], 10.871940554611)]
check("PiT(10,Pe_c) = -8.40e-15",
      float(crit.loc[crit.k.idxmax(), "Pi_T"]), -8.40e-15, 1e-16)
true("no positive sample at Pe=10.9", bool((sub.Pi_T < 0).all()))
roots = pd.read_csv(os.path.join(DENSE, "ultimate_roots.csv"))
r85 = roots[np.isclose(roots.Pe, 10.85)].iloc[0]
check("interior root k* = 9.338218553141", float(r85.k_root), 9.338218553141,
      1e-9, rel=True)
check("|PiT(k*)| = 1.92e-20", float(r85.abs_Pi_root), 1.92e-20, 1e-22)
true("root type at Pe=10.85 is interior", r85.root_type, "interior_root")
true("root type at Pe=Pe_c is upper boundary",
     roots[np.isclose(roots.Pe, 10.871940554611)].root_type.iloc[0],
     "upper_boundary_root")
loc = pd.read_csv(os.path.join(DENSE, "ultimate_local_boundary_check.csv"))
check("PiT(10,Pe_c-1e-4) = +5.00641e-9",
      float(loc[np.isclose(loc.delta_Pe, -1e-4)].Pi_T.iloc[0]), 5.00641e-9, 1e-13)
check("PiT(10,Pe_c+1e-4) = -5.00648e-9",
      float(loc[np.isclose(loc.delta_Pe, 1e-4)].Pi_T.iloc[0]), -5.00648e-9, 1e-13)
check("endpoint slope -5.01e-5",
      (float(loc[np.isclose(loc.delta_Pe, 1e-4)].Pi_T.iloc[0])
       - float(loc[np.isclose(loc.delta_Pe, -1e-4)].Pi_T.iloc[0])) / 2e-4,
      -5.0064e-5, 1e-8)

# --------------------------------------------------- DNS estimator identity (B7)
DNSDIR = os.path.join(ROOT, "dns_physical_validation_suite")
sym = pd.read_csv(os.path.join(DNSDIR, "dns_symmetry_test.csv")).set_index("case")
base = sym.loc["baseline_asymmetric"]
check("eps_a*eps_s*2 = 3.0e-4", 2 * 0.03 * 0.005, 3.0e-4, 0)
check("PiT_DNS = (q_+-q_-)/(2 eps_a eps_s)",
      base.q_difference / (2 * 0.03 * 0.005), base.PiT, 1e-12, rel=True)
check("q_+-q_- = -1.40289979e-5", base.q_difference, -1.40289979e-5, 1e-12)
check("reflected q_+-q_- = +1.40289979e-5",
      sym.loc["reflected_asymmetric"].q_difference, +1.40289979e-5, 1e-12)
check("reflection flux residual 3.18e-17",
      abs(sym.loc["reflected_asymmetric"].q_difference
          + base.q_difference), 3.18e-17, 1e-18)
check("reflection PiT residual 1.06e-13",
      abs(sym.loc["reflected_asymmetric"].PiT - base.PiT), 1.06e-13, 1e-14)
check("symmetric-forcing q_+-q_- = 8.13e-18",
      sym.loc["symmetric"].q_difference, 8.13e-18, 1e-19)
check("wall-flux conservation q_++q_- = -1.9e-7",
      base.q_sum, -1.9e-7, 1e-8)
check("PiT_DNS baseline = -4.6763e-2", base.PiT, -4.6763e-2, 1e-5, rel=True)
mas = pd.read_csv(os.path.join(DNSDIR, "dns_master_validation.csv"))
re5 = mas[mas.test.str.contains("Re=5")].iloc[0]
check("A_DNS(Re=5) = -1.5878e-3", re5.value, -1.5878e-3, 1e-7)

# ------------------------------------------------ merged Re_c locus (7 points)
loc7 = pd.read_csv(os.path.join(HERE, "D1_Re_c_all.csv"))
true("D1_Re_c_all resolves 7 crossings", loc7.shape[0], 7)
for k, want in [(0.4, 79.080214), (0.6, 51.252439)]:
    check(f"Re_c(k={k}) low-k extension",
          float(loc7.loc[loc7.k == k, "Re_c"].iloc[0]), want, 5e-7)
lowk = pd.read_csv(os.path.join(HERE, "D1_Re_grid_lowk.csv"))
true("low-k grid rows", lowk.shape[0], 20)
true("low-k grid reaches Re=120", float(lowk.Re.max()), 120.0)
true("low-k antisymmetry residual < 1e-17", float(lowk.antisym_err.max()) < 1e-17)
fit7 = json.load(open(os.path.join(HERE, "D1_Re_c_fit.json")))
check("seven-point fit prefactor 21.50", fit7["C"], 21.50, 0.01)
check("seven-point fit exponent -1.47", fit7["n"], -1.47, 0.01)
check("seven-point fit max deviation 11.0 %",
      100 * fit7["max_rel_deviation"], 11.0, 0.1)
check("five-point fit max deviation 7.7 %",
      100 * fit7["five_point"]["max_rel_deviation"], 7.729, 0.02)

# ------------------------------------------ DNS deviation model of the table
sweep = pd.read_csv(os.path.join(DNSDIR, "dns_re_sweep.csv"))
pert = {(round(r.k, 6), round(r.Re, 6)): r.Pi_T for r in grid.itertuples()}
flat, rel = [], []
for _, r in sweep.iterrows():
    p = pert[(0.8, round(float(r.Re), 6))]
    e = float(r.PiT) - p
    if abs(p) < 1e-3:
        flat.append(abs(e))
    elif abs(p) > 1e-2:
        rel.append(abs(e / p))
check("flat DNS deviation 2.71e-4", max(flat), 2.71e-4, 1e-6)
check("relative DNS deviation 1.52 %", max(rel), 1.52e-2, 1e-5)

# ------------------------------------------------ DNS pending-prediction table
pred = pd.read_csv(os.path.join(HERE, "DNS_D345_predictions.csv"))
true("prediction table has 13 cases", pred.shape[0], 13)
for cid, want in [("D3_k0.8_Re20", -1.928365e-6), ("D3_k0.8_Re40", 3.242275e-7),
                  ("D3_k1.0_Re15", -1.828239e-6), ("D3_k1.0_Re25", 4.671537e-7),
                  ("D3_k1.2_Re10", -3.702864e-6), ("D3_k1.2_Re20", 7.719747e-7),
                  ("D3b_k0.6_Re40", -3.886322e-7), ("D3b_k0.6_Re70", 1.740707e-7),
                  ("D3c_k0.8_Re30_ea02", -2.350614e-7),
                  ("D4_k0.8_Re50", 3.989744e-7),
                  ("D5_phi_plus_pi_over_2", 8.005465e-5),
                  ("D5_phi_minus_pi_over_2", -8.005465e-5)]:
    got = float(pred.loc[pred.case_id == cid, "q_diff_predicted"].iloc[0])
    check(f"predicted q_+-q_- {cid}", got, want, 5e-12)
check("D5 phi=0 predicts exact cancellation",
      abs(float(pred.loc[pred.case_id == "D5_phi_0",
                         "q_diff_predicted"].iloc[0])), 0.0, 1e-30)
for r in pred.itertuples():
    rec = 2.0 * (r.c_bilinear * r.PiT_perturbative
                 + r.d_bilinear * r.PiT_perp_perturbative)
    check(f"phase-law identity {r.case_id}", rec, r.q_diff_predicted, 1e-16)
    if not np.isfinite(r.margin_H1):
        continue          # phi=0 carries no bilinear amplitude: no resolution
    dev = max(2.71e-4, 1.52e-2 * abs(r.PiT_projected))
    check(f"resolution H1 {r.case_id}", abs(r.PiT_projected) / dev,
          r.margin_H1, 1e-9)
strong = pred[(pred.family != "D.3c") & pred.margin_H1.notna()]
true("all non-ladder cases sign-resolvable (H1 >= 2)",
     bool((strong.margin_H1 >= 2.0).all()))
ladder = pred[pred.family == "D.3c"].iloc[0]
true("ladder case discriminates the models",
     bool(ladder.margin_H1 < 2.0 <= ladder.margin_H2))

# ------------------------------------------------ DNS campaign measurements
mdf = pd.read_csv(os.path.join(HERE, "DNS_D345_measured.csv"))
m = mdf.set_index("case_id")
QUOTED_Q = {
    "D3_k0.8_Re20": -1.799329e-6, "D3_k0.8_Re40": +3.809729e-7,
    "D3_k1.0_Re15": -1.604171e-6, "D3_k1.0_Re25": +5.910363e-7,
    "D3_k1.2_Re10": -3.326895e-6, "D3_k1.2_Re20": +9.454173e-7,
    "D3b_k0.6_Re40": -3.542744e-7, "D3b_k0.6_Re70": +1.910334e-7,
    "D3c_k0.8_Re30_ea02": +9.007819e-8, "D4_k0.8_Re50": +4.416458e-7,
    "D5_phi_plus_pi_over_2": +8.115924e-5,
    "D5_phi_minus_pi_over_2": -8.115932e-5}
true("DNS campaign has 13 runs", len(mdf), 13)
for cid, want in QUOTED_Q.items():
    check("measured q_diff " + cid, m.loc[cid, "q_diff_dns"], want, 5e-7,
          rel=True)
check("measured q_diff D5_phi_0 < 1e-17", abs(m.loc["D5_phi_0", "q_diff_dns"]),
      0.0, 1e-17)

QUOTED_DEV = {
    "D3_k0.8_Re20": 4.30e-4, "D3_k0.8_Re40": 1.89e-4,
    "D3_k1.0_Re15": 7.47e-4, "D3_k1.0_Re25": 4.13e-4,
    "D3_k1.2_Re10": 1.25e-3, "D3_k1.2_Re20": 5.78e-4,
    "D3b_k0.6_Re40": 1.15e-4, "D3b_k0.6_Re70": 5.66e-5,
    "D3c_k0.8_Re30_ea02": 2.71e-4, "D4_k0.8_Re50": 1.42e-4,
    "D5_phi_plus_pi_over_2": 1.23e-3, "D5_phi_minus_pi_over_2": 1.23e-3}
for cid, want in QUOTED_DEV.items():
    check("coefficient deviation " + cid, m.loc[cid, "deviation"], want,
          5e-3, rel=True)
QUOTED_RATIO = {"D3_k0.8_Re20": 1.6, "D3_k1.0_Re15": 2.8,
                "D3_k1.2_Re10": 4.6, "D3_k1.0_Re25": 1.5,
                "D3_k1.2_Re20": 2.1}
for cid, want in QUOTED_RATIO.items():
    check("delta_Pi multiple " + cid, m.loc[cid, "dev_over_model"], want, 0.06)
bound = m[m.dev_over_model.notna()]
check("cases within pre-registered bound",
      int((bound.dev_over_model <= 1.0).sum()), 7, 0)
check("cases exceeding bound", int((bound.dev_over_model > 1.0).sum()), 5, 0)

side = m[m.family.isin(["D.3", "D.3b"])]
check("bracket runs at four wavenumbers", len(side), 8, 0)
true("every bracket run returns the predicted sign",
     bool((np.sign(side.q_diff_dns) == np.sign(side.q_diff_predicted)).all()))
true("brackets cover k = 0.6, 0.8, 1.0, 1.2",
     sorted(side.k.unique().tolist()), [0.6, 0.8, 1.0, 1.2])
nzd = m[m.q_diff_predicted.abs() > 1e-12]
check("sign agreement among non-degenerate runs",
      int((np.sign(nzd.q_diff_dns) == np.sign(nzd.q_diff_predicted)).sum()),
      11, 0)
check("total campaign wall time (min)", mdf.wall_seconds.sum() / 60.0, 24.8,
      0.05)
check("largest final residual", float(mdf.final_residual.max()), 2.2e-8,
      2.2e-8)
true("per-case runtimes inside 90-144 s",
     bool(mdf.wall_seconds.min() > 90.0 and mdf.wall_seconds.max() < 145.0))

lad = m.loc["D3c_k0.8_Re30_ea02"]
check("ladder deviation over absolute floor", lad.deviation / 2.71e-4, 1.00,
      0.01)
true("ladder sign disagrees with prediction",
     bool(np.sign(lad.q_diff_dns) != np.sign(lad.q_diff_predicted)))
d4 = m.loc["D4_k0.8_Re50"]
check("D4 measured PiT", d4.PiT_dns, 1.472e-3, 1e-3, rel=True)
check("D4 predicted PiT", d4.PiT_pred, 1.330e-3, 1e-3, rel=True)
check("D4 deviation over floor", d4.deviation / 2.71e-4, 0.52, 0.01)

qp, qm = m.loc["D5_phi_plus_pi_over_2"], m.loc["D5_phi_minus_pi_over_2"]
check("quadrature Nu_20 measured", qp.q_diff_dns / 2.0, 4.0580e-5, 1e-4,
      rel=True)
check("quadrature relative deviation", qp.deviation / abs(qp.PiT_pred), 1.38e-2,
      1e-3, rel=True)
check("quadrature antisymmetry residual", abs(qp.q_diff_dns + qm.q_diff_dns),
      8.01e-11, 1e-12)
check("quadrature antisymmetry relative",
      abs(qp.q_diff_dns + qm.q_diff_dns) / abs(qp.q_diff_dns), 9.9e-7, 1e-8)
true("quadrature runs have opposite sign",
     bool(np.sign(qp.q_diff_dns) != np.sign(qm.q_diff_dns)))

per_k = bound[bound.family.isin(["D.3", "D.3b"])].groupby("k").deviation.max()
true("deviation grows monotonically with k",
     bool(per_k.is_monotonic_increasing))
for kk, want in [(0.6, 1.15e-4), (0.8, 4.30e-4), (1.0, 7.47e-4),
                 (1.2, 1.25e-3)]:
    check(f"max deviation at k={kk}", float(per_k.loc[kk]), want, 5e-3,
          rel=True)

ana = json.load(open(os.path.join(HERE, "DNS_D345_analysis.json")))
check("analysis reversal agreements", ana["reversal"]["n_sign_agreements"], 8, 0)
check("analysis within-bound count", ana["deviation_model"]["n_within_bound"],
      7, 0)
check("analysis worst ratio", ana["deviation_model"]["max_dev_over_model"],
      4.624, 0.002)
true("analysis worst case is k=1.2",
     ana["deviation_model"]["worst_case"] == "D3_k1.2_Re10")
true("analysis ladder reading is H1", ana["ladder"]["reading"] == "H1")
check("analysis total wall minutes", ana["total_wall_minutes"], 24.8, 0.05)
true("multi-wavenumber DNS figure exists",
     os.path.exists(os.path.join(HERE, "FIG_v2_10_dns_multik.png")))

# ------------------------------------------- units, parity and k-scaling (B3)
udf = pd.read_csv(os.path.join(HERE, "DNS_D345_measured.csv"))
conv = {cid: 2.0 * float(udf[udf.case_id == cid].amp_bilinear.iloc[0])
        for cid in ["D3_k0.8_Re20", "D3c_k0.8_Re30_ea02",
                    "D5_phi_plus_pi_over_2"]}
check("conversion factor, reference amplitude", conv["D3_k0.8_Re20"],
      3.0e-4, 1e-12)
check("conversion factor, ladder amplitude", conv["D3c_k0.8_Re30_ea02"],
      1.2e-3, 1e-12)
check("conversion factor, equal-amplitude quadrature",
      conv["D5_phi_plus_pi_over_2"], 9.0e-4, 1e-12)
dq = {cid: float(udf[udf.case_id == cid].q_dev.iloc[0])
      for cid in udf.case_id}
check("ladder flux gap 3.25e-7", dq["D3c_k0.8_Re30_ea02"], 3.25e-7, 5e-9)
check("reference-amplitude flux gap 1.29e-7 at k=0.8, Re=20",
      dq["D3_k0.8_Re20"], 1.29e-7, 5e-9)
check("ladder flux gap divided by its factor is the quoted deviation",
      dq["D3c_k0.8_Re30_ea02"] / conv["D3c_k0.8_Re30_ea02"], 2.709e-4,
      5e-3, rel=True)
par = udf[udf.case_id == "D5_phi_0"].iloc[0]
check("parity flux difference 3.1e-18", float(par.q_diff_dns), 3.1e-18, 5e-20)
check("parity wall flux -1.145796985e-7", float(par.q_top), -1.145796985e-7,
      1e-9, rel=True)
check("parity relative cancellation 2.7e-11",
      abs(par.q_diff_dns) / abs(par.q_top), 2.74e-11, 5e-13)
true("parity run has vanishing antisymmetric amplitude",
     bool(par.eps_a == 0.0))
brk = udf[udf.family.isin(["D.3", "D.3b"])].dropna(subset=["deviation"])
true("k-scaling fit uses the eight bracket runs", brk.shape[0], 8)
xk, yk = np.log10(brk.k.values), np.log10(brk.deviation.values)
slope, inter = np.polyfit(xk, yk, 1)
resk = yk - (slope * xk + inter)
sig = np.sqrt((resk ** 2).sum() / (len(xk) - 2))
err_slope = sig * np.sqrt(len(xk) / ((xk - xk.mean()) ** 2).sum())
r2k = 1 - (resk ** 2).sum() / ((yk - yk.mean()) ** 2).sum()
check("remainder exponent k^3.4", slope, 3.411, 0.01)
check("remainder exponent uncertainty 1.7", err_slope, 1.72, 0.02)
check("remainder fit R^2 0.84", r2k, 0.840, 0.005)
true("deviation grows monotonically with k (maxima)",
     bool(np.all(np.diff(brk.groupby("k").deviation.max().values) > 0)))
sc = pd.read_csv(os.path.join(ROOT, "dns_physical_validation_suite",
                              "dns_spatial_convergence.csv"))
true("convergence ladder spans 1776 to 7008 degrees of freedom",
     sorted(int(v) for v in sc.resolution), [1776, 3136, 4880, 7008])
check("spatial drift at the production mesh 3.1e-4",
      float(sc[sc.resolution == 3136].abs_error_vs_finest.iloc[0]), 3.11e-4,
      5e-2, rel=True)
check("production mesh degrees of freedom is Nx*Ny", 64 * 49, 3136, 0)
true("campaign residuals below the quoted 2.2e-8",
     bool(udf.final_residual.max() < 2.2e-8))

# ------------------------------------------------ DNS Grashof pair (Cor. 1 B2)
gdf = pd.read_csv(os.path.join(HERE, "DNS_D2_grpair_measured.csv"))
gj = json.load(open(os.path.join(HERE, "DNS_D2_grpair_analysis.json")))
true("Gr-pair DNS matrix has 4 runs", gdf.shape[0], 4)
true("Gr-pair at k=0.8 and Pe=7",
     bool((gdf.k == 0.8).all() and (gdf.Pe == 7.0).all()))
true("Gr-pair uses the Grashof pair 50 and 200",
     sorted(gdf.Gr.unique().tolist()), [50.0, 200.0])
true("Gr-pair straddles Re_c = 30.838", sorted(gdf.Re.unique().tolist()),
     [20.0, 50.0])
for cid, want in [("D2_Gr50_Re20", -8.996526e-07),
                  ("D2_Gr200_Re20", -3.598839e-06),
                  ("D2_Gr50_Re50", +2.207512e-07),
                  ("D2_Gr200_Re50", +8.833327e-07)]:
    check("Gr-pair %s measured q_+-q_-" % cid,
          float(gdf[gdf.case_id == cid].q_diff_dns.iloc[0]), want, 5e-7,
          rel=True)
gp = {p["Re"]: p for p in gj["pairs"]}
check("measured Gr-ratio at Re=20 is 4.0003", gp[20.0]["Gr_ratio_measured"],
      4.0003, 5e-5)
check("measured Gr-ratio at Re=50 is 4.0015", gp[50.0]["Gr_ratio_measured"],
      4.0015, 5e-5)
check("worst Gr-ratio relative error 3.7e-4",
      gj["worst_Gr_ratio_rel_error"], 3.71e-4, 2e-6)
true("Gr-pair signs agree within each pair", gj["all_signs_agree"])
true("Gr-pair negative below Re_c", gp[20.0]["sign"] == "negative")
true("Gr-pair positive above Re_c", gp[50.0]["sign"] == "positive")
for Re in (20.0, 50.0):
    check("deviation ratio at Re=%.0f is proportional to Gr" % Re,
          gp[Re]["dev_ratio"], 4.0, 0.02)
    sub = gdf[gdf.Re == Re]
    check("A relative error at Re=%.0f is Gr-independent" % Re,
          sub.A_rel_err.max() / sub.A_rel_err.min() - 1.0, 0.0, 5e-3)
check("relative error of A at Re=20 is 6.69%",
      gdf[gdf.Re == 20].A_rel_err.mean(), 6.69e-2, 5e-3, rel=True)
check("relative error of A at Re=50 is 10.7%",
      gdf[gdf.Re == 50].A_rel_err.mean(), 0.107, 5e-3, rel=True)
check("Gr-pair total wall time 491 s", gj["total_wall_seconds"], 491.0, 2.0)
true("Gr-pair residuals below 1.0e-8",
     bool(gj["max_final_residual"] < 1.0e-8))
true("Gr-pair runtimes within 100-150 s",
     bool(100.0 <= gj["runtime_range"][0] and gj["runtime_range"][1] <= 150.0))
true("all Gr-pair margins pass the gate", bool((gdf.margin >= 2.0).all()))
check("weakest retained margin is 2.5", gdf.margin.min(), 2.5, 0.05)
for Re, Gr, wm, wp in [(20.0, 50.0, 11.9, 0.3500), (20.0, 200.0, 47.4, 0.3500),
                       (50.0, 50.0, 2.5, 0.1400), (50.0, 200.0, 9.8, 0.1400)]:
    row = gdf[(gdf.Re == Re) & (gdf.Gr == Gr)].iloc[0]
    check("tab:Gr_dns margin at Re=%.0f, Gr=%.0f" % (Re, Gr), row.margin,
          wm, 0.05)
    check("tab:Gr_dns Pr at Re=%.0f, Gr=%.0f" % (Re, Gr), row.Pr, wp, 5e-4)
true("Richardson numbers quoted in the caveat are correct",
     bool(abs(200.0 / 20.0 ** 2 - 0.50) < 1e-12
          and abs(200.0 / 50.0 ** 2 - 0.08) < 1e-12))
rej = {r["Re"]: r["weakest"] for r in gj["rejected"]}
check("rejected Re=40 design has margin 1.99", rej[40.0], 1.99, 0.01)
check("rejected Re=32 design has margin 0.43", rej[32.0], 0.43, 0.01)

# ------------------------------------------------- Gr-pair figure (Cor. 1)
fig_path = os.path.join(HERE, "FIG_v2_11_dns_grpair.png")
true("Gr-pair figure exists", os.path.exists(fig_path))
fj = json.loads(open(os.path.join(HERE, "FIG_v2_11_dns_grpair.json"),
                     encoding="utf-8").read())
true("perturbative Gr branches collapse below 2e-10",
     bool(fj["branch_collapse_rel"] < 2.0e-10))
check("figure crossing line is the perturbative Re_c", fj["Re_c"],
      30.8384965, 5e-6)
true("figure grid has 25 points per branch", fj["grid_points"] == 25)
true("figure marker offset is the one in the caption",
     fj["marker_offset_re"] == 0.7)
check("figure pair agreement in A at Re=20", fj["pair_agreement_rel"]["20"],
      6.342e-5, 2e-8)
check("figure pair agreement in A at Re=50", fj["pair_agreement_rel"]["50"],
      3.711e-4, 5e-8)
for Re, Gr, aval in [(20.0, 50.0, -5.997684e-05), (20.0, 200.0, -5.998064e-05),
                     (50.0, 50.0, 1.471675e-05), (50.0, 200.0, 1.472221e-05)]:
    row = gdf[(gdf.Re == Re) & (gdf.Gr == Gr)].iloc[0]
    check("fig:dns_grpair A_dns at Re=%.0f Gr=%.0f" % (Re, Gr), row.A_dns,
          aval, 5e-11)

# ------------------------------------------------------- supplementary package
man = pd.read_csv(os.path.join(ROOT, "supplementary_material", "MANIFEST.csv"))
true("supplementary MANIFEST has 64 files", man.shape[0], 64)
true("every manifest file exists",
     all(os.path.exists(os.path.join(ROOT, "supplementary_material",
                                     g["group"], f))
         for g, f in zip(man.to_dict("records"), man.file)))

print(f"\n{len(fails)} failure(s)" + (": " + "; ".join(fails) if fails else ""))
sys.exit(1 if fails else 0)

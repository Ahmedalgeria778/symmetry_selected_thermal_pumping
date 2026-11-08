#!/usr/bin/env python3
"""Unwrapped phase-lag analysis of the Re-driven reversal (plan item C.2).

For each Re we record the continuous phase difference between the antisymmetric
velocity mode and the symmetric thermal mode, and locate the phase crossing that
coincides with the Pi_T reversal.
"""
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = Path(r"C:\Users\user\Documents\Documents - Copie\mixed_conv_verification")
OUT = WORK / "campaign_v2_symmetry_selection"
_spec = importlib.util.spec_from_file_location("v4", WORK / "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py")
v4 = importlib.util.module_from_spec(_spec); sys.modules["v4"] = v4; _spec.loader.exec_module(v4)

K, PE, GR = 0.8, 7.0, 100.0
RES = np.array([5., 6., 7., 8., 9., 10., 12., 15., 18., 20., 22., 25., 27., 28., 29.,
                30., 30.5, 30.838497, 31., 32., 33., 35., 38., 40., 45., 50.])
Y = np.array([0.0, 0.25, 0.5, 0.75])

rows = []
t0 = time.time()
for Re in RES:
    ts = v4.solve_thermal_bvp(K, PE, True)
    ta = v4.solve_thermal_bvp(K, PE, False)
    hs = v4.solve_hydro_bvp(K, Re, GR, ts)
    ha = v4.solve_hydro_bvp(K, Re, GR, ta)
    y = np.linspace(-1., 1., 1001)
    Ts, Ta = ts.theta(y), ta.theta(y)
    Ts_y, Ta_y = ts.theta_y(y), ta.theta_y(y)
    Us, Ua, Vs, Va = hs.U(y), ha.U(y), hs.V(y), ha.V(y)
    Z1 = Ua * np.conjugate(-1j * K * Ts) + Va * np.conjugate(Ts_y)
    Z2 = Us * np.conjugate(-1j * K * Ta) + Vs * np.conjugate(Ta_y)
    X = 0.5 * (np.real(Z1) + np.real(Z2))
    m = v4.solve_mean_bvp(y, X, PE)
    thy = m.theta_y(y)
    pi_t = 0.5 * (-thy[-1] - thy[0])
    # individual term contributions
    terms = {"T_Ua_dXs": 0.5*np.real(Ua*np.conjugate(-1j*K*Ts)),
             "T_Us_dXa": 0.5*np.real(Us*np.conjugate(-1j*K*Ta)),
             "T_Va_dYs": 0.5*np.real(Va*np.conjugate(Ts_y)),
             "T_Vs_dYa": 0.5*np.real(Vs*np.conjugate(Ta_y))}
    contrib = {}
    for nm, prof in terms.items():
        mm = v4.solve_mean_bvp(y, prof, PE)
        tt = mm.theta_y(y)
        contrib[nm] = 0.5 * (-tt[-1] - tt[0])
    # phases at selected wall-normal stations
    ths, tha, vs, va = ts.theta(Y), ta.theta(Y), hs.V(Y), ha.V(Y)
    dphi_s = np.angle(vs / ths)
    dphi_a = np.angle(va / tha)
    # overlap factors: cos of the local phase difference entering each kernel term
    ov_a_dXs = np.real(Ua * np.conjugate(-1j*K*Ts)) / (np.abs(Ua)*np.abs(K*Ts) + 1e-300)
    ov_a_dYs = np.real(Va * np.conjugate(Ts_y)) / (np.abs(Va)*np.abs(Ts_y) + 1e-300)
    rows.append(dict(Re=Re, Pi_T=pi_t,
                     **contrib,
                     dphi_s_0=dphi_s[0], dphi_s_5=dphi_s[2],
                     dphi_a_0=dphi_a[0], dphi_a_5=dphi_a[2],
                     abs_Th_s_0=np.abs(ths[0]), abs_V_s_0=np.abs(vs[0]),
                     abs_Th_a_0=np.abs(tha[0]), abs_V_a_0=np.abs(va[0]),
                     ov_a_dXs_0=float(ov_a_dXs[len(y)//2]),
                     ov_a_dYs_0=float(ov_a_dYs[len(y)//2]),
                     ratio_a_dYs=float(np.max(np.abs(Va*np.conjugate(Ts_y))) /
                                       (np.max(np.abs(Ua*np.conjugate(-1j*K*Ts))) + 1e-300))))
    print("Re=%6.2f Pi_T=% .6e  T_Ua=% .4e T_Va=% .4e  dphi_s0=%+.4f  dphi_a0=%+.4f"
          % (Re, pi_t, contrib["T_Ua_dXs"], contrib["T_Va_dYs"], dphi_s[0], dphi_a[0]), flush=True)

df = pd.DataFrame(rows)
# unwrap the centreline phases along the Re sequence
for col in ("dphi_s_0", "dphi_s_5", "dphi_a_0", "dphi_a_5"):
    df[col + "_unw"] = np.unwrap(df[col].values)
df["sum_terms"] = df[["T_Ua_dXs", "T_Us_dXa", "T_Va_dYs", "T_Vs_dYa"]].sum(axis=1)
df["term_balance_err"] = df.sum_terms - df.Pi_T
df.to_csv(OUT / "C2_phase_unwrapped.csv", index=False)
print("max |sum(terms) - Pi_T| = %.3e" % df.term_balance_err.abs().max(), flush=True)

fig, ax = plt.subplots(2, 1, figsize=(6.4, 7.2), sharex=True)
a = ax[0]
a.plot(df.Re, df.T_Ua_dXs, "o-", label=r"$\mathcal{L}[\langle U^{(a)}\partial_X\Theta^{(s)}\rangle]$")
a.plot(df.Re, df.T_Va_dYs, "s-", label=r"$\mathcal{L}[\langle V^{(a)}\partial_Y\Theta^{(s)}\rangle]$")
a.plot(df.Re, df.T_Us_dXa, "^-", ms=4, label=r"$\mathcal{L}[\langle U^{(s)}\partial_X\Theta^{(a)}\rangle]$")
a.plot(df.Re, df.T_Vs_dYa, "v-", ms=4, label=r"$\mathcal{L}[\langle V^{(s)}\partial_Y\Theta^{(a)}\rangle]$")
a.axhline(0, color="k", lw=0.8)
a.axvline(30.838497, color="grey", ls="--", lw=0.9)
a.set_ylabel("contribution to $\Pi_T$")
a.legend(fontsize=7)
a.set_title(r"$k=0.8$, $Pe=7$, $Gr=100$: term-wise balance and reversal at $Re_c=30.838$")
b = ax[1]
b.plot(df.Re, df.dphi_s_0_unw, "o-", label=r"$\Delta\varphi^{(s)}=\arg(\widehat V^{(s)}/\Theta^{(s)})$ at $Y=0$")
b.plot(df.Re, df.dphi_a_0_unw, "s-", label=r"$\Delta\varphi^{(a)}=\arg(\widehat V^{(a)}/\Theta^{(a)})$ at $Y=0$")
b.axvline(30.838497, color="grey", ls="--", lw=0.9)
b.axhline(0, color="k", lw=0.8)
b.set_xlabel(r"$Re$ at fixed $Pe=7$")
b.set_ylabel(r"unwrapped phase lag [rad]")
b.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "FIG_v2_08_phase_crossing.png", dpi=200)
plt.close(fig)

# where does each phase lag cross zero / pi/2?
def cross(col):
    x = df.Re.values; yv = df[col].values
    out = []
    for i in range(len(x)-1):
        if np.sign(yv[i]) != np.sign(yv[i+1]) and np.isfinite(yv[i]) and np.isfinite(yv[i+1]):
            t = yv[i]/(yv[i]-yv[i+1])
            out.append(float(x[i] + t*(x[i+1]-x[i])))
    return out

df["T_Va_dYs_plus"] = df.T_Ua_dXs + df.T_Va_dYs
rep = dict(Re_c_Pi_T=30.838497,
           zero_cross_dphi_s=cross("dphi_s_0_unw"),
           zero_cross_dphi_a=cross("dphi_a_0_unw"),
           term_T_Va_dYs_plus_T_Ua_zero=cross("T_Va_dYs_plus"),
           max_balance_err=float(df.term_balance_err.abs().max()))
df["T_Va_dYs_plus"] = df.T_Ua_dXs + df.T_Va_dYs
rep["term_balance_Ua_Va_zero"] = cross("T_Va_dYs_plus")
json.dump(rep, open(OUT / "C2_phase_crossings.json", "w"), indent=2, default=float)
print(json.dumps(rep, indent=2, default=float), flush=True)
print("TOTAL %.0f s" % (time.time()-t0))

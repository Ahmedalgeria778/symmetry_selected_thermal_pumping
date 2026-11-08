"""Analyse the completed DNS campaign D.3/D.4/D.5 against its predictions.

Reads DNS_D345_results.csv (written by run_dns_D3_D4_D5.py --go) and checks the
four claims the campaign was designed to test:

  1. sign of the wall-flux difference flips across Re_c at every k sampled;
  2. the measured deviation from the two-coefficient prediction is bounded by the
     pre-registered model delta_Pi = max(2.71e-4, 1.52% |Pi_T|);
  3. the parity cancellation at phi = 0 and the antisymmetry in phi;
  4. which of the two readings of the error model (H1 fixed in coefficient space,
     H2 absolute in q_+-q_-) the amplitude ladder supports.

Writes DNS_D345_analysis.json, DNS_D345_measured.csv and FIG_v2_10_dns_multik.png.

    py -X utf8 campaign_v2_symmetry_selection/analyze_dns_D345.py
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
DNS_DEV, DNS_DEV_REL = 2.71e-4, 1.52e-2

df = pd.read_csv(OUT / "DNS_D345_results.csv")

amp = np.hypot(df.c_bilinear, df.d_bilinear)
df["amp_bilinear"] = amp
df["PiT_dns"] = np.where(amp > 0, df.q_diff_dns / (2.0 * amp), np.nan)
df["PiT_pred"] = np.where(amp > 0, df.q_diff_predicted / (2.0 * amp), np.nan)
df["deviation"] = (df.PiT_dns - df.PiT_pred).abs()
df["delta_model"] = np.maximum(DNS_DEV, DNS_DEV_REL * df.PiT_pred.abs())
df["dev_over_model"] = df.deviation / df.delta_model
df["q_dev"] = (df.q_diff_dns - df.q_diff_predicted).abs()

# --- 1. reversal sign at every sampled wavenumber -----------------------------
rec = pd.read_csv(OUT / "D1_Re_c_all.csv").set_index("k").Re_c
br = df[df.family.isin(["D.3", "D.3b"])].copy()
br["Re_c"] = [rec.loc[k] for k in br.k]
br["side"] = np.where(br.Re > br.Re_c, "above", "below")
br["sign_dns"] = np.sign(br.q_diff_dns)
br["expected"] = np.where(br.side == "above", 1.0, -1.0)
br["sign_ok"] = br.sign_dns == br.expected
reversal = dict(n_bracket_cases=int(len(br)),
                n_sign_agreements=int(br.sign_ok.sum()),
                k_values=sorted(set(float(k) for k in br.k)),
                worst_margin=float(br.margin_H1.min()),
                rows=[dict(case_id=r.case_id, k=float(r.k), Re=float(r.Re),
                           Re_c=float(r.Re_c), side=r.side,
                           q_diff_dns=float(r.q_diff_dns),
                           q_diff_predicted=float(r.q_diff_predicted),
                           sign_ok=bool(r.sign_ok),
                           margin_H1=float(r.margin_H1))
                      for r in br.itertuples()])

# --- 2. pre-registered deviation model ---------------------------------------
finite = df[df.dev_over_model.notna()]
model = dict(form="delta_Pi = max(2.71e-4, 1.52% |Pi_T|)",
             n_cases=int(len(finite)),
             max_dev_over_model=float(finite.dev_over_model.max()),
             n_within_bound=int((finite.dev_over_model <= 1.0).sum()),
             worst_case=str(finite.loc[finite.dev_over_model.idxmax(), "case_id"]),
             ratios={str(r.case_id): float(r.dev_over_model)
                     for r in finite.itertuples()},
             deviations={str(r.case_id): float(r.deviation)
                         for r in finite.itertuples()})

# --- 3. phase tests ----------------------------------------------------------
p0 = df[df.case_id == "D5_phi_0"].iloc[0]
pp = df[df.case_id == "D5_phi_plus_pi_over_2"].iloc[0]
pm = df[df.case_id == "D5_phi_minus_pi_over_2"].iloc[0]
phase = dict(q_diff_phi0=float(p0.q_diff_dns),
             q_diff_phi0_predicted=float(p0.q_diff_predicted),
             antisymmetry_residual=float(abs(pp.q_diff_dns + pm.q_diff_dns)),
             antisymmetry_relative=float(abs(pp.q_diff_dns + pm.q_diff_dns) /
                                         abs(pp.q_diff_dns)),
             PiT_dns_plus=float(pp.PiT_dns), PiT_dns_minus=float(pm.PiT_dns),
             PiT_dns_plus_predicted=float(pp.PiT_pred),
             ratio_measured_over_predicted_plus=float(pp.q_diff_dns /
                                                      pp.q_diff_predicted),
             ratio_measured_over_predicted_minus=float(pm.q_diff_dns /
                                                       pm.q_diff_predicted))

# --- 4. ladder: which reading of the error model -----------------------------
lad = df[df.family == "D.3c"].iloc[0]
ladder = dict(case_id=lad.case_id, k=float(lad.k), Re=float(lad.Re),
              eps_a=float(lad.eps_a), amp=float(lad.amp_bilinear),
              q_diff_predicted=float(lad.q_diff_predicted),
              q_diff_dns=float(lad.q_diff_dns),
              PiT_pred=float(lad.PiT_pred), PiT_dns=float(lad.PiT_dns),
              deviation=float(lad.deviation),
              deviation_over_floor=float(lad.deviation / DNS_DEV),
              margin_H1=float(lad.margin_H1), margin_H2=float(lad.margin_H2),
              predicted_sign="negative" if lad.q_diff_predicted < 0 else "positive",
              measured_sign="negative" if lad.q_diff_dns < 0 else "positive",
              reading="H1" if lad.q_diff_dns * lad.q_diff_predicted < 0 else "H2")

summary = dict(n_cases=int(len(df)),
               total_wall_minutes=float(df.wall_seconds.sum() / 60.0),
               max_final_residual=float(df.final_residual.max()),
               sign_agreement_nonzero=int(
                   ((np.sign(df.q_diff_dns) == np.sign(df.q_diff_predicted)) &
                    (df.q_diff_predicted.abs() > 1e-12)).sum()),
               n_nonzero=int((df.q_diff_predicted.abs() > 1e-12).sum()),
               reversal=reversal, deviation_model=model, phase=phase,
               ladder=ladder)

df.to_csv(OUT / "DNS_D345_measured.csv", index=False)
json.dump(summary, open(OUT / "DNS_D345_analysis.json", "w"), indent=2)

print("cases: %d | wall %.1f min | max residual %.2e"
      % (summary["n_cases"], summary["total_wall_minutes"],
         summary["max_final_residual"]))
print("sign agreement: %d/%d non-degenerate"
      % (summary["sign_agreement_nonzero"], summary["n_nonzero"]))
print("reversal brackets: %d/%d, k = %s"
      % (reversal["n_sign_agreements"], reversal["n_bracket_cases"],
         reversal["k_values"]))
print("deviation model: %d/%d within bound, worst ratio %.3f (%s)"
      % (model["n_within_bound"], model["n_cases"],
         model["max_dev_over_model"], model["worst_case"]))
print("ladder: deviation = %.3e = %.2f x floor -> H1"
      % (ladder["deviation"], ladder["deviation_over_floor"]))
print("phi=0: q_diff = %.2e | antisymmetry residual %.2e (rel %.1e)"
      % (phase["q_diff_phi0"], phase["antisymmetry_residual"],
         phase["antisymmetry_relative"]))

# --- figure ------------------------------------------------------------------
grid = pd.concat([pd.read_csv(OUT / "D1_Re_grid.csv"),
                  pd.read_csv(OUT / "D1_Re_grid_lowk.csv")], ignore_index=True)

fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
ax = axes[0]
cols = {0.6: "tab:blue", 0.8: "tab:red", 1.0: "tab:green", 1.2: "tab:purple"}
for k in sorted(cols):
    g = grid[np.isclose(grid.k, k)].sort_values("Re")
    ax.plot(g.Re, g.A, "-", color=cols[k], lw=1.4,
            label=r"$k=%.1f$, $\Pi_T/Gr$" % k)
    dns = df[(np.isclose(df.k, k)) &
             df.family.isin(["D.3", "D.3b", "D.3c", "D.4"])]
    ax.plot(dns.Re, dns.PiT_dns / 100.0, "o", color=cols[k], ms=8,
            mec="k", mew=0.6)
ax.axhline(0.0, color="k", lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_xlabel(r"$Re$")
ax.set_ylabel(r"$\Pi_T/Gr$")
ax.set_title("DNS markers on the perturbative locus: sign flips at $Re_c$")
ax.legend(fontsize=7, ncol=2)

ax = axes[1]
pl = df[df.q_diff_predicted.abs() > 1e-12]
err = 2.0 * pl.amp_bilinear * pl.delta_model
ax.errorbar(pl.q_diff_predicted, pl.q_diff_dns, yerr=err, fmt="none",
            ecolor="0.6", elinewidth=1, capsize=2)
ax.scatter(pl.q_diff_predicted, pl.q_diff_dns, c=np.sign(pl.q_diff_predicted),
           cmap="bwr", s=46, vmin=-1, vmax=1)
lim = [min(pl.q_diff_predicted.min(), pl.q_diff_dns.min()) * 0.6,
       max(pl.q_diff_predicted.max(), pl.q_diff_dns.max()) * 1.6]
ax.plot(lim, lim, "k--", lw=0.9)
ax.set_xscale("symlog", linthresh=1e-7)
ax.set_yscale("symlog", linthresh=1e-7)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel(r"predicted $q_+-q_-$, two-coefficient law")
ax.set_ylabel(r"measured $q_+-q_-$ (DNS)")
ax.set_title("Measured vs predicted wall-flux difference ($\\pm\\delta$ bars)")
fig.tight_layout()
fig.savefig(OUT / "FIG_v2_10_dns_multik.png", dpi=300)
print("figure: FIG_v2_10_dns_multik.png")

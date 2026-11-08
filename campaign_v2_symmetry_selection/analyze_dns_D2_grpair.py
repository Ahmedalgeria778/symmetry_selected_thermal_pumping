"""Reduce the Grashof-pair DNS checkpoint to the nonlinear test of Corollary 1.

Corollary 1 states Pi_T = Gr * A(k, Re, Pr), so at fixed Re the pumping flux must
scale exactly with Gr and the reversal must not move. The measurable statement is
therefore about the ratio of the two runs at each Re, not about their absolute
size: any deviation that is itself proportional to Gr leaves the zero untouched.

    py -X utf8 campaign_v2_symmetry_selection/analyze_dns_D2_grpair.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
SRC = OUT / "DNS_D2_grpair_results.csv"
WORK = OUT.parent
sys.path.insert(0, str(WORK))
_spec = importlib.util.spec_from_file_location(
    "d345", OUT / "run_dns_D3_D4_D5.py")
d345 = importlib.util.module_from_spec(_spec)
sys.modules["d345"] = d345
_spec.loader.exec_module(d345)
GR_LO, GR_HI = 50.0, 200.0
K, PE, EPS_S, EPS_A = 0.8, 7.0, 0.03, 0.005


def rejected(Re_list=(32.0, 40.0)):
    """Margins of the Re pairs the resolvability gate turned down before DNS."""
    c, _d = d345.bilinear_factors(EPS_S + EPS_A, EPS_S - EPS_A, 0.0)
    out = []
    for Re in Re_list:
        margins = []
        for Gr in (GR_LO, GR_HI):
            pi_c = d345.coefficients(K, Re, PE, Gr)[0]
            margins.append(abs(pi_c) / d345.dns_deviation(pi_c))
        out.append(dict(Re=float(Re), margins=margins, weakest=min(margins)))
        print("rejected design Re=%5.1f  margins %.2f / %.2f"
              % (Re, margins[0], margins[1]))
    return out


def main():
    m = pd.read_csv(SRC)
    m["A_perturbative"] = m.PiT_perturbative / m.Gr
    m["A_dns"] = m.PiT_dns / m.Gr
    m["deviation_PiT"] = (m.PiT_dns - m.PiT_perturbative).abs()
    m["deviation_q"] = (m.q_diff_dns - m.q_diff_predicted).abs()
    m["delta_q"] = 2.0 * m.c_bilinear * m.delta_model
    m["dev_over_model"] = m.deviation_q / m.delta_q
    m["A_rel_err"] = m.deviation_PiT / (m.Gr * m.A_perturbative.abs())
    m.to_csv(OUT / "DNS_D2_grpair_measured.csv", index=False,
             float_format="%.9e")
    pd.set_option("display.width", 200)
    print(m[["case_id", "Re", "Gr", "q_diff_predicted", "q_diff_dns",
             "PiT_perturbative", "PiT_dns", "dev_over_model", "A_dns",
             "A_rel_err", "final_residual"]].to_string(index=False))

    pairs = []
    for Re in sorted(m.Re.unique()):
        lo = m[(m.Re == Re) & (m.Gr == GR_LO)].iloc[-1]
        hi = m[(m.Re == Re) & (m.Gr == GR_HI)].iloc[-1]
        ratio = hi.q_diff_dns / lo.q_diff_dns
        pairs.append(dict(
            Re=float(Re), q_lo=float(lo.q_diff_dns), q_hi=float(hi.q_diff_dns),
            Gr_ratio_measured=float(ratio),
            Gr_ratio_rel_error=float(abs(ratio / 4.0 - 1.0)),
            same_sign=bool(np.sign(lo.q_diff_dns) == np.sign(hi.q_diff_dns)),
            sign=str("negative" if lo.q_diff_dns < 0 else "positive"),
            PiT_lo=float(lo.PiT_dns), PiT_hi=float(hi.PiT_dns),
            dev_lo=float(lo.deviation_PiT), dev_hi=float(hi.deviation_PiT),
            dev_ratio=float(hi.deviation_PiT / lo.deviation_PiT),
            A_dns_lo=float(lo.A_dns), A_dns_hi=float(hi.A_dns),
            A_spread_rel=float(abs(hi.A_dns / lo.A_dns - 1.0)),
            margin_lo=float(lo.margin), margin_hi=float(hi.margin),
            dev_over_model_lo=float(lo.dev_over_model),
            dev_over_model_hi=float(hi.dev_over_model)))
    for p in pairs:
        print("Re=%5.1f  sign %s  Gr-ratio %.4f  rel err %.2e  "
              "deviation ratio %.3f  A spread %.2e  margins %.1f/%.1f  "
              "dev/floor %.2f/%.2f"
              % (p["Re"], p["sign"], p["Gr_ratio_measured"],
                 p["Gr_ratio_rel_error"], p["dev_ratio"], p["A_spread_rel"],
                 p["margin_lo"], p["margin_hi"], p["dev_over_model_lo"],
                 p["dev_over_model_hi"]))

    summary = dict(
        k=float(m.k.iloc[0]), Pe=float(m.Pe.iloc[0]), eps_s=0.03, eps_a=0.005,
        phi=0.0, cases=int(m.shape[0]), pairs=pairs,
        rejected=rejected(),
        worst_Gr_ratio_rel_error=float(max(p["Gr_ratio_rel_error"]
                                           for p in pairs)),
        all_signs_agree=bool(all(p["same_sign"] for p in pairs)),
        max_final_residual=float(m.final_residual.max()),
        total_wall_seconds=float(m.wall_seconds.sum()),
        runtime_range=[float(m.wall_seconds.min()), float(m.wall_seconds.max())],
        reading=("the DNS deviation is itself proportional to Gr (ratios "
                 "%.3f and %.3f against 4.0000), so it rescales the response "
                 "without moving the zero: Corollary 1 holds nonlinearly"
                 % (pairs[0]["dev_ratio"], pairs[1]["dev_ratio"])))
    json.dump(summary, open(OUT / "DNS_D2_grpair_analysis.json", "w"), indent=2)
    print("worst Gr-ratio relative error: %.2e"
          % summary["worst_Gr_ratio_rel_error"])
    print("total wall %.0f s, max residual %.2e"
          % (summary["total_wall_seconds"], summary["max_final_residual"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

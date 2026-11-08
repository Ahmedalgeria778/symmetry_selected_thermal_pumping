"""Nonlinear check of Corollary 1: the reversal point does not move with Gr.

Corollary 1 says Pi_T = Gr * A(k, Re, Pr), so at fixed Re the pumping coefficient
is exactly proportional to Gr and the crossing Re_c is Gr-independent. The
perturbative verification already gives a spread of 1.2e-10 in Re_c over five
Grashof numbers. This script extends the check to the nonlinear regime with a
2x2 DNS matrix: Re on both sides of Re_c(0.8) = 30.838, and Gr = 50, 200.

Two things must hold. At each Re the measured ratio
q_+-q_-(200) / q_+-q_-(50) must be 4 (linearity in Gr, DNS), and the sign must
be negative at Re = 20 and positive at Re = 50 for both Gr (the crossing has not
moved with Gr). The dry run rejects any design whose weakest run falls inside
the DNS deviation measured in the D.3/D.4/D.5 campaign; Re = 40 was rejected by
that gate because its Gr = 50 member has a margin of only 1.99.

    py -X utf8 campaign_v2_symmetry_selection/run_dns_D2_grpair.py
    py -X utf8 campaign_v2_symmetry_selection/run_dns_D2_grpair.py --go
"""
import importlib.util
import os
import sys
import time
from pathlib import Path

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

K, PE, EPS_S, EPS_A = 0.8, 7.0, 0.03, 0.005
GR_LO, GR_HI = 50.0, 200.0
RE_LO, RE_HI = 20.0, 50.0
CHECKPOINT = OUT / "DNS_D2_grpair_results.csv"


def predict(Re, Gr):
    c, _d = d345.bilinear_factors(EPS_S + EPS_A, EPS_S - EPS_A, 0.0)
    pi_c, _pi_p, anti_c, _anti_p = d345.coefficients(K, Re, PE, Gr)
    delta = d345.dns_deviation(pi_c)
    return dict(Re=Re, Gr=Gr, Pr=PE / Re, c=c, PiT=pi_c,
                q_diff_predicted=2.0 * c * pi_c, delta_model=delta,
                margin=abs(pi_c) / delta, antisym_err=anti_c)


def main():
    go = "--go" in sys.argv
    Re_list, Gr_list = [RE_LO, RE_HI], [GR_LO, GR_HI]
    for a in sys.argv[1:]:
        if a.startswith("--Re="):
            Re_list = [float(v) for v in a.split("=")[1].split(",")]

    rows = [predict(Re, Gr) for Re in Re_list for Gr in Gr_list]
    print("k=%.1f  Pe=%.1f  Pr=Pe/Re  eps_s=%.3f eps_a=%.3f  phi=0"
          % (K, PE, EPS_S, EPS_A))
    for r in rows:
        print("  Re=%5.1f Gr=%6.1f  PiT=%+.6e  q_+-q_-=%+.6e  delta=%.2e"
              "  margin=%5.1f  antisym=%.1e"
              % (r["Re"], r["Gr"], r["PiT"], r["q_diff_predicted"],
                 r["delta_model"], r["margin"], r["antisym_err"]))
    for Re in Re_list:
        pair = [r for r in rows if r["Re"] == Re]
        print("  predicted Gr-ratio at Re=%.1f: %.6f"
              % (Re, pair[1]["q_diff_predicted"] / pair[0]["q_diff_predicted"]))
    worst = min(r["margin"] for r in rows)
    if worst < 2.0:
        print("NOT RESOLVABLE: the weakest run sits inside the DNS deviation"
              " measured at k=0.8. Move the Re pair further from Re_c=30.838"
              " before integrating.")
        return 1
    if not go:
        print("DRY RUN: nothing integrated. Re-invoke with --go.")
        return 0

    from dns_physical_validation_suite import run_dns  # noqa: E402

    try:
        done = set(pd.read_csv(CHECKPOINT).case_id)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        done = set()
    handle = open(CHECKPOINT, "a", newline="")
    new_file = not done
    for r in rows:
        cid = "D2_Gr%.0f_Re%.0f" % (r["Gr"], r["Re"])
        if cid in done:
            print("skip", cid)
            continue
        t0 = time.time()
        res = run_dns(eps1=EPS_S + EPS_A, eps2=EPS_S - EPS_A, k=K, Re=r["Re"],
                      Pr=r["Pr"], Gr=r["Gr"], Nx=64, Ny=49, dt=0.0015,
                      steps=30000, phi=0.0, verbose=False)
        q = res["q_top"] - res["q_bottom"]
        pd.DataFrame([dict(case_id=cid, k=K, Re=r["Re"], Pr=r["Pr"], Pe=PE,
                           Gr=r["Gr"], eps1=EPS_S + EPS_A,
                           eps2=EPS_S - EPS_A, c_bilinear=r["c"],
                           PiT_perturbative=r["PiT"],
                           q_diff_predicted=r["q_diff_predicted"],
                           delta_model=r["delta_model"], margin=r["margin"],
                           q_top=res["q_top"], q_bottom=res["q_bottom"],
                           q_diff_dns=q, PiT_dns=q / (2.0 * r["c"]),
                           final_residual=res["final_residual"],
                           wall_seconds=time.time() - t0)]).to_csv(
            handle, mode="a", header=new_file, index=False)
        new_file = False
        handle.flush()
        os.fsync(handle.fileno())
        print("DONE %-16s q_diff=%+.6e (predicted %+.6e) resid=%.2e %.0fs"
              % (cid, q, r["q_diff_predicted"], res["final_residual"],
                 time.time() - t0), flush=True)
    handle.close()

    m = pd.read_csv(CHECKPOINT)
    for Re in sorted(m.Re.unique()):
        lo = m[(m.Re == Re) & (m.Gr == GR_LO)].iloc[-1]
        hi = m[(m.Re == Re) & (m.Gr == GR_HI)].iloc[-1]
        print("Re=%.1f measured Gr-ratio %.4f (Corollary 1: 4.0000),"
              " same sign: %s"
              % (Re, hi.q_diff_dns / lo.q_diff_dns,
                 bool(np.sign(lo.q_diff_dns) == np.sign(hi.q_diff_dns))))
    return 0


if __name__ == "__main__":
    sys.exit(main())

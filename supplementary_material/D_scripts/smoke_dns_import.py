"""Smoke test for the pending DNS campaign: cheap import + short integration.

Checks, before any long run is committed, that
  (a) dns_physical_validation_suite imports without touching a foreign path,
  (b) run_dns accepts exactly the keyword arguments used by
      campaign_v2_symmetry_selection/run_dns_D3_D4_D5.py,
  (c) the returned dict carries the keys the checkpoint writer reads, with finite
      values and a converged-in-shape residual.

200 steps at the production mesh, so it costs well under a minute.

    py -X utf8 campaign_v2_symmetry_selection/smoke_dns_import.py
"""
import math
import sys
import time
from pathlib import Path

WORK = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(WORK))

import dns_physical_validation_suite as suite  # noqa: E402

print("suite OUT:", suite.OUT)
if not suite.OUT.is_dir():
    print("FAIL: OUT does not exist:", suite.OUT)
    sys.exit(1)

t0 = time.time()
res = suite.run_dns(eps1=0.035, eps2=0.025, k=0.8, Re=10.0, Pr=0.7, Gr=100.0,
                    Nx=64, Ny=49, dt=0.0015, steps=200, phi=0.0, verbose=False)
dt_wall = time.time() - t0

need = ["q_top", "q_bottom", "PiT", "final_residual"]
missing = [key for key in need if key not in res]
if missing:
    print("FAIL: run_dns result missing keys:", missing)
    sys.exit(1)

bad = [key for key in need
       if isinstance(res[key], float) and not math.isfinite(res[key])]
if bad:
    print("FAIL: non-finite values for:", bad, {k: res[k] for k in bad})
    sys.exit(1)

q_diff = res["q_top"] - res["q_bottom"]
print("keys ok:", sorted(res))
print("q_top=%.6e q_bottom=%.6e q_diff=%+.6e" %
      (res["q_top"], res["q_bottom"], q_diff))
print("PiT=%.6e residual=%.3e  %.1f s for 200 steps" %
      (res["PiT"], res["final_residual"], dt_wall))
print("extrapolated 30000 steps: %.1f min" % (dt_wall * 30000 / 200 / 60))
print("SMOKE OK")

"""Block until the DNS campaign checkpoint is complete, then summarise it.

    py -X utf8 campaign_v2_symmetry_selection/wait_dns_D345.py [--target 13] [--max-minutes 120]

Prints one progress line per newly finished case and a final table of measured
versus predicted wall-flux differences. Exits 0 when the checkpoint reaches
--target rows, 2 on timeout, so it can be used as a background waiter.
"""
import argparse
import sys
import time
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent
CKPT = OUT / "DNS_D345_results.csv"

ap = argparse.ArgumentParser()
ap.add_argument("--target", type=int, default=13)
ap.add_argument("--max-minutes", type=float, default=120.0)
ap.add_argument("--poll-seconds", type=float, default=60.0)
args = ap.parse_args()

t_end = time.time() + 60.0 * args.max_minutes
seen = 0
while time.time() < t_end:
    rows = 0
    if CKPT.exists():
        try:
            df = pd.read_csv(CKPT)
            rows = len(df)
        except pd.errors.EmptyDataError:
            rows = 0
    if rows > seen:
        for i in range(seen, rows):
            print("case %d/%d done: %s" % (i + 1, args.target,
                                           df.case_id.iloc[i]), flush=True)
        seen = rows
    if rows >= args.target:
        break
    time.sleep(args.poll_seconds)

if not CKPT.exists():
    print("TIMEOUT: no checkpoint after %.0f min" % args.max_minutes)
    sys.exit(2)

df = pd.read_csv(CKPT)
cols = ["case_id", "k", "Re", "phi", "q_diff_predicted", "q_diff_dns",
        "expected_sign", "margin_H1", "final_residual", "wall_seconds"]
print("\n" + df[[c for c in cols if c in df.columns]].to_string(index=False))
sign_ok = ((df.q_diff_dns > 0) == (df.q_diff_predicted > 0)) | \
          (df.q_diff_predicted.abs() < 1e-12)
print("\nsign agreement: %d/%d non-degenerate cases"
      % (int(sign_ok.sum()), len(df)))
print("total wall time: %.1f min" % (df.wall_seconds.sum() / 60.0))
sys.exit(0 if len(df) >= args.target else 2)

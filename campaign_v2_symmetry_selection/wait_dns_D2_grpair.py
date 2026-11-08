"""Block until the Gr-pair DNS matrix reports its summary, then print the logs."""
import time
from pathlib import Path

OUT = Path(__file__).resolve().parent
LOG = OUT / "log_DNS_D2_grpair.txt"
ERR = OUT / "log_DNS_D2_grpair.err.txt"
MARKER = "measured Gr-ratio"
DEADLINE = time.time() + 2400
STALL = 600

last = -1
while time.time() < DEADLINE:
    text = LOG.read_text(encoding="utf-8", errors="replace") if LOG.exists() else ""
    size = LOG.stat().st_size if LOG.exists() else 0
    if size != last:
        print("log %d bytes, DONE lines %d" % (size, text.count("DONE")),
              flush=True)
        last = size
    if MARKER in text:
        break
    if last and time.time() - LOG.stat().st_mtime > STALL:
        print("STALLED: no log write for %d s" % STALL, flush=True)
        break
    time.sleep(45)

print("=" * 70)
print(text)
if ERR.exists() and ERR.stat().st_size:
    print("--- STDERR ---")
    print(ERR.read_text(encoding="utf-8", errors="replace")[-2000:])

"""Reduced Pe_c endpoint figure (plan item D.6).

Panel (a): the response Pi_T(k, Pe) for the two retained off-endpoint
calculations plus the endpoint itself. Panel (b): zoom on the high-wavenumber
region where the zero-crossing branch migrates to the upper boundary. Only
data from the CPU-authoritative dense scan already present in this project are
used; the three post-endpoint Pe values of the original scan are not plotted.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Pe_c = 10.871940554611
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "validation_outputs_ULTIMATE_DENSE_CPU",
                   "ultimate_dense_cpu_scan.csv")
OUT = os.path.join(HERE, "FIG_v2_09_Pe_c_endpoint.png")

KEPT = [10.85, Pe_c, 10.9]
LABELS = {10.85: r"$Pe=10.85$", Pe_c: r"$Pe=Pe_c=10.8719406$",
          10.9: r"$Pe=10.9$"}
COLORS = {10.85: "#1f77b4", Pe_c: "#444444", 10.9: "#d62728"}

df = pd.read_csv(SRC)
curves = {pe: df[np.isclose(df["Pe"], pe)].sort_values("k") for pe in KEPT}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=200)

for pe in KEPT:
    c = curves[pe]
    ax1.plot(c["k"], c["Pi_T"], color=COLORS[pe], lw=1.8, label=LABELS[pe])
ax1.axhline(0.0, color="k", lw=0.8, ls="--")
ax1.plot([9.338218553141], [0.0], "o", ms=6, color="#1f77b4")
ax1.plot([10.0], [0.0], "s", ms=6, color="#444444")
ax1.set_xlabel(r"$k$")
ax1.set_ylabel(r"$\Pi_T(k,Pe)$")
ax1.set_title(r"(a) full admissible interval, $Re=10$, $Gr=100$")
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

for pe in KEPT:
    c = curves[pe]
    m = c["k"] >= 8.0
    ax2.plot(c.loc[m, "k"], c.loc[m, "Pi_T"], color=COLORS[pe], lw=2.0,
             label=LABELS[pe])
ax2.axhline(0.0, color="k", lw=0.8, ls="--")
hi = max(curves[pe].loc[curves[pe]["k"] >= 8.0, "Pi_T"].max() for pe in KEPT)
lo = min(curves[pe].loc[curves[pe]["k"] >= 8.0, "Pi_T"].min() for pe in KEPT)
ax2.set_xlim(8.0, 10.05)
ax2.set_ylim(1.15 * lo, 1.15 * hi)
ax2.set_xlabel(r"$k$")
ax2.set_ylabel(r"$\Pi_T(k,Pe)$")
ax2.set_title(r"(b) zoom: interior zero $\to$ boundary contact $\to$ no zero")
ax2.annotate(r"$k^*=9.3382185531$", xy=(9.3382, 0.0),
             xytext=(8.15, lo + 0.18 * (hi - lo)), fontsize=8,
             arrowprops=dict(arrowstyle="->", lw=0.9))
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT)
print("wrote", OUT)
for pe in KEPT:
    c = curves[pe]
    print(f"Pe={pe}: min={c['Pi_T'].min():.6e} max={c['Pi_T'].max():.6e}")

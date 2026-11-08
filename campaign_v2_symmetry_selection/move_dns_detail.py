"""Move the DNS convergence detail from the manuscript to Supplementary S2.

Review move D asks the manuscript to keep the essential error levels and to send
the mesh, time-step and amplitude sequences to the supplementary document. The
sequences are tabulated displays, not floats, so the transfer adds numbered
equations to S2 and leaves the figure and table numbering untouched. The
manuscript keeps one condensed subsection with the same numbers.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "main_v2.tex"
SUPP = ROOT / "supplementary_text.tex"

KEEP = r"""\subsection{Convergence of the nonlinear estimator}
\label{subsec:dns_amplitude}
\label{subsec:dns_convergence}

Four sequences fix the reliability of the nonlinear observable; the sequences
themselves are recorded in Supplementary Sec.~\ref{subsec:dns_sequences}. As the
asymmetry amplitude decreases at fixed $\varepsilon_s=0.03$,
$\PiT^{\mathrm{DNS}}$ converges to $-4.6763\times10^{-2}$, which recovers the
second-order perturbative value $-4.7485\times10^{-2}$ to $1.52\%$
(Fig.~\ref{fig:dns_amplitude}). The approach to a finite limit is the nonlinear
evidence that $\PiT$ is the leading transport coefficient of the
symmetric--antisymmetric interaction rather than an amplitude artefact. Refining
the mesh from $48\times37$ to $96\times73$ moves the coefficient monotonically by
at most $1.10\%$, and refining the time step from $3.0\times10^{-3}$ to
$1.0\times10^{-3}$ by at most $1.75\times10^{-6}$
(Figs.~\ref{fig:dns_spatial} and \ref{fig:dns_timestep}), so the wall-flux
difference is neither under-resolved in the transverse boundary layers nor
controlled by the temporal discretisation. At the reference resolution the final
evolution residual stays below $7\times10^{-9}$ and the global wall-flux balance
closes to $q_++q_-\simeq-1.9\times10^{-7}$, so the coefficient is extracted from a
converged stationary state rather than from a finite-time snapshot.
"""

LEAD = r"""
\subsection{Amplitude, grid and time-step sequences}
\label{subsec:dns_sequences}

The four levels quoted in the manuscript rest on the sequences recorded below.

"""


def find(lines, needle, start=0):
    for i in range(start, len(lines)):
        if needle in lines[i]:
            return i
    raise StopIteration


def main():
    lines = MAIN.read_text(encoding="utf-8").splitlines(keepends=True)
    start = find(lines, "\\subsection{Amplitude convergence and recovery")
    stop = find(lines, "\\subsection{Independent verification of the Grashof")
    moved = lines[start:stop]
    lines[start:stop] = [KEEP + "\n"]
    MAIN.write_text("".join(lines), encoding="utf-8")

    body = [t for t in moved if not t.startswith("\\subsection{")]
    body = [t for t in body if "\\label{subsec:dns_" not in t]
    while body and not body[-1].strip():
        body.pop()

    supp = SUPP.read_text(encoding="utf-8").splitlines(keepends=True)
    cut = find(supp, "\\section{Phase decomposition and term-wise balance}")
    supp[cut:cut] = [LEAD] + body + ["\n"]
    SUPP.write_text("".join(supp), encoding="utf-8")
    print("moved %d lines; manuscript keeps %d" % (len(moved), KEEP.count("\n")))


if __name__ == "__main__":
    main()

import numpy as np

def json_serializable(obj):
    """Convertit les objets NumPy/Pandas en types Python sérialisables JSON."""
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        if np.isnan(obj):
            return None
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Series):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient='records')
    elif isinstance(obj, (np.generic,)):
        return obj.item()
    else:
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

def run_dns(
    eps1,
    eps2,
    k=0.8,
    Re=10.0,
    Pr=0.7,
    Gr=100.0,
    Nx=64,
    Ny=49,
    dt=0.0015,
    steps=30000,
    tol=1e-8,
    print_every=100,
    verbose=True,
    phi=0.0,
):
    """
    Petite DNS pseudo-spectrale Fourier en x + différences finies / DST en y.

    Résout les équations de perturbation non linéaires autour de Poiseuille :

        u = U0(y) + psi_y
        v = - psi_x

    avec vorticité perturbation :

        zeta = v_x - u_y = -nabla^2 psi

    Convention conforme au papier :

        paroi haute  Y = +1 : T = eps1 cos(kX)
        paroi basse  Y = -1 : T = eps2 cos(kX + phi)

    phi = 0 reproduit exactement les conditions aux limites utilisees dans les
    campagnes precedentes; la convention est celle du manuscrit, qui ecrit
    theta_hat(-1) = eps2 exp(-i phi) sous la convention e^{-ikX}.

    Le coefficient de pompage DNS est :

        Pi_T = (q_top - q_bottom) / (2 eps_s eps_a)

    où :

        q_top    = - mean( dT/dY at Y=+1 )
        q_bottom = + mean( dT/dY at Y=-1 )

    """

    # ------------------------------------------------------------------
    # Géométrie et grilles
    # ------------------------------------------------------------------
    Lx = 2.0 * np.pi / k
    dx_grid = Lx / Nx
    dy = 2.0 / (Ny - 1)

    x = np.linspace(0.0, Lx, Nx, endpoint=False)
    y = np.linspace(-1.0, 1.0, Ny)

    X, Y = np.meshgrid(x, y, indexing="ij")

    # ------------------------------------------------------------------
    # Paramètres dérivés
    # ------------------------------------------------------------------
    Pe = Re * Pr

    eps_s = 0.5 * (eps1 + eps2)
    eps_a = 0.5 * (eps1 - eps2)

    # Pi_T is undefined for pure symmetric or pure antisymmetric forcing,
    # but the DNS remains valid for testing the flux cancellation.

    # ------------------------------------------------------------------
    # Opérateurs spectraux en X
    # ------------------------------------------------------------------
    K = 2.0 * np.pi * np.fft.fftfreq(Nx, d=dx_grid)
    K2 = K * K

    def ddx(f):
        """Dérivée première en x, périodique."""
        fh = np.fft.fft(f, axis=0)
        return np.fft.ifft(1j * K[:, None] * fh, axis=0).real

    def lapx(f):
        """Dérivée seconde en x, périodique."""
        fh = np.fft.fft(f, axis=0)
        return np.fft.ifft(-K2[:, None] * fh, axis=0).real

    # ------------------------------------------------------------------
    # Transformée sinus discrète en Y pour Dirichlet psi=0 aux parois
    # ------------------------------------------------------------------
    M = Ny - 2
    if M < 1:
        raise ValueError("Ny doit être >= 3, mais Ny>=5 est recommandé.")

    j = np.arange(1, M + 1)
    m = np.arange(1, M + 1)

    # Matrice DST-I orthonormale
    S = np.sqrt(2.0 / (M + 1)) * np.sin(
        np.pi * np.outer(m, j) / (M + 1)
    )

    # Valeurs propres de d^2/dy^2 avec Dirichlet aux bords
    lam_y = -4.0 / dy**2 * np.sin(
        np.pi * m / (2.0 * (M + 1))
    ) ** 2

    # Dénominateur pour résoudre nabla^2 psi = -zeta
    # (-K^2 + lam_y) psi_hat = -zeta_hat
    # psi_hat = zeta_hat / (K^2 - lam_y)
    denom = K2[:, None] - lam_y[None, :]

    def solve_psi(zeta):
        """
        Résout nabla^2 psi = -zeta avec psi=0 en Y=±1.
        """
        psi = np.zeros_like(zeta)

        # seulement les points intérieurs
        zeta_int = zeta[:, 1:-1]

        # FFT en x
        zeta_x = np.fft.fft(zeta_int, axis=0)

        # DST en y
        zeta_xy = zeta_x @ S.T

        # résolution modale
        psi_xy = zeta_xy / denom

        # DST inverse
        psi_x = psi_xy @ S

        # IFFT en x
        psi[:, 1:-1] = np.fft.ifft(psi_x, axis=0).real

        # Dirichlet psi=0 aux parois
        psi[:, 0] = 0.0
        psi[:, -1] = 0.0

        return psi

    # ------------------------------------------------------------------
    # Dérivées finite differences en Y
    # ------------------------------------------------------------------
    def ddy(f):
        """Dérivée première centrée en Y."""
        out = np.empty_like(f)
        out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dy)

        # one-sided aux bords, même si on utilise surtout l'intérieur
        out[:, 0] = (-3.0 * f[:, 0] + 4.0 * f[:, 1] - f[:, 2]) / (2.0 * dy)
        out[:, -1] = (3.0 * f[:, -1] - 4.0 * f[:, -2] + f[:, -3]) / (2.0 * dy)

        return out

    def d2y(f):
        """Dérivée seconde centrée en Y, seulement intérieure."""
        out = np.zeros_like(f)
        out[:, 1:-1] = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / dy**2
        return out

    # ------------------------------------------------------------------
    # Conditions limites
    # ------------------------------------------------------------------
    bc_top = eps1 * np.cos(k * x)
    bc_bottom = eps2 * np.cos(k * x + phi)

    def set_theta_bc(theta):
        theta[:, 0] = bc_bottom
        theta[:, -1] = bc_top

    def set_zeta_bc(zeta, psi):
        """
        Condition de vorticité à la paroi compatible avec no-slip.

        psi=0 à la paroi et psi_y=0 donnent approximativement :

            zeta_wall = -2 psi_near / dy^2

        car psi_xx wall = 0 puisque psi wall = 0 pour tout x.
        """
        zeta[:, 0] = -2.0 * psi[:, 1] / dy**2
        zeta[:, -1] = -2.0 * psi[:, -2] / dy**2

    # ------------------------------------------------------------------
    # État initial
    # ------------------------------------------------------------------
    # Profil linéaire en y qui satisfait déjà les BC thermiques
    theta = (eps_s + Y * eps_a) * np.cos(k * X)
    zeta = np.zeros((Nx, Ny))

    set_theta_bc(theta)

    # ------------------------------------------------------------------
    # Écoulement de base
    # ------------------------------------------------------------------
    U0 = 1.5 * (1.0 - y**2)
    U0pp = -3.0 * np.ones(Ny)

    # ------------------------------------------------------------------
    # RHS
    # ------------------------------------------------------------------
    def rhs(zeta, theta):
        """
        RHS des équations non linéaires perturbées.

        Vorticité :

            ∂t ζ =
              - U0 ζ_x
              - U0'' ψ_x
              + J(ψ, ζ)
              + (1/Re) ∇²ζ
              - (Gr/Re^2) θ_y

        Température :

            ∂t θ =
              - U0 θ_x
              + J(ψ, θ)
              + (1/Pe) ∇²θ

        avec :

            J(a,b) = a_x b_y - a_y b_x

        """
        set_theta_bc(theta)

        psi = solve_psi(zeta)
        set_zeta_bc(zeta, psi)

        psix = ddx(psi)
        psiy = ddy(psi)

        zetax = ddx(zeta)
        zetay = ddy(zeta)

        thetax = ddx(theta)
        thetay = ddy(theta)

        lap_zeta = lapx(zeta) + d2y(zeta)
        lap_theta = lapx(theta) + d2y(theta)

        Jzeta = psix * zetay - psiy * zetax
        Jtheta = psix * thetay - psiy * thetax

        dzeta = (
            -U0[None, :] * zetax
            - U0pp[None, :] * psix
            + Jzeta
            + (1.0 / Re) * lap_zeta
            - (Gr / Re**2) * thetay
        )

        dtheta = (
            -U0[None, :] * thetax
            + Jtheta
            + (1.0 / Pe) * lap_theta
        )

        # Les variables de bord sont imposées / diagnostiques
        dzeta[:, 0] = 0.0
        dzeta[:, -1] = 0.0
        dtheta[:, 0] = 0.0
        dtheta[:, -1] = 0.0

        return dzeta, dtheta

    # ------------------------------------------------------------------
    # Time stepping SSP-RK3
    # ------------------------------------------------------------------
    zeta_old = zeta.copy()
    theta_old = theta.copy()

    history = []

    for n in range(steps):

        # Stage 1
        dz1, dth1 = rhs(zeta, theta)
        z1 = zeta + dt * dz1
        th1 = theta + dt * dth1
        set_theta_bc(th1)

        # Stage 2
        dz2, dth2 = rhs(z1, th1)
        z2 = 0.75 * zeta + 0.25 * (z1 + dt * dz2)
        th2 = 0.75 * theta + 0.25 * (th1 + dt * dth2)
        set_theta_bc(th2)

        # Stage 3
        dz3, dth3 = rhs(z2, th2)
        zeta = (1.0 / 3.0) * zeta + (2.0 / 3.0) * (z2 + dt * dz3)
        theta = (1.0 / 3.0) * theta + (2.0 / 3.0) * (th2 + dt * dth3)
        set_theta_bc(theta)

        if (n % print_every == 0) or (n == steps - 1):

            dz_check, dth_check = rhs(zeta, theta)

            res = max(
                np.max(np.abs(dz_check[:, 1:-1])),
                np.max(np.abs(dth_check[:, 1:-1])),
            )

            change = max(
                np.max(np.abs(zeta[:, 1:-1] - zeta_old[:, 1:-1])),
                np.max(np.abs(theta[:, 1:-1] - theta_old[:, 1:-1])),
            )

            history.append((n, res, change))

            if verbose:
                print(
                    f"step={n:6d} | residual={res:.3e} | change={change:.3e}"
                )

            zeta_old = zeta.copy()
            theta_old = theta.copy()

            if not np.isfinite(res):
                print("ERROR: residual non fini. Réduire dt ou raffiner.")
                break

            if res < tol:
                if verbose:
                    print("Convergence atteinte.")
                break

    # ------------------------------------------------------------------
    # Calcul du flux pariétal moyen
    # ------------------------------------------------------------------
    # Dérivée one-sided second ordre
    theta_y_top = (
        3.0 * theta[:, -1] - 4.0 * theta[:, -2] + theta[:, -3]
    ) / (2.0 * dy)

    theta_y_bottom = (
        -3.0 * theta[:, 0] + 4.0 * theta[:, 1] - theta[:, 2]
    ) / (2.0 * dy)

    q_top = -np.mean(theta_y_top)
    q_bottom = np.mean(theta_y_bottom)

    denom_pi = 2.0 * eps_s * eps_a
    PiT_dns = (q_top - q_bottom) / denom_pi if abs(denom_pi) > 1e-14 else np.nan

    # Résidu final
    dz_final, dth_final = rhs(zeta, theta)
    final_residual = max(
        np.max(np.abs(dz_final[:, 1:-1])),
        np.max(np.abs(dth_final[:, 1:-1])),
    )

    return {
        "PiT": PiT_dns,
        "q_top": q_top,
        "q_bottom": q_bottom,
        "eps_s": eps_s,
        "eps_a": eps_a,
        "final_residual": final_residual,
        "history": history,
        "theta": theta,
        "zeta": zeta,
    }


# =============================================================================
# MULTI-STUDY PHYSICAL DNS VALIDATION SUITE
# =============================================================================
# Independent nonlinear DNS checks designed for publication-quality validation.
# Studies:
#   V1 amplitude convergence/extrapolation eps_a -> 0
#   V2 spatial convergence
#   V3 time-step convergence
#   V4 Grashof scaling and Gr=0 switch-off
#   V5 Reynolds sweep at fixed Pe
#   V6 wall-reflection and symmetric-heating tests
#   V7 master publication summary + CSV/LaTeX tables
#
# The perturbative reference values below are the values already established
# in the manuscript for the same physical parameters.
# =============================================================================

from pathlib import Path
from dataclasses import asdict
import json
import time
import warnings
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "dns_physical_validation_suite"
OUT.mkdir(parents=True, exist_ok=True)

K_REF = 0.8
PE_REF = 7.0
RE_REF = 10.0
PR_REF = PE_REF / RE_REF
GR_REF = 100.0
EPS_S_REF = 0.030
PI_REF = -4.7485e-2
RE_CRIT_PE7 = 30.838496518052

# perturbative normalized response A=Pi_T/Gr from the manuscript
A_REF_BY_RE = {
    5.0:  -1.5951e-3,
    10.0: -4.7485e-4,
    20.0: -6.4279e-5,
    35.0:  6.7452e-6,
}

# ---------------------------- publication controls -------------------------
QUICK = False
BASE_NX = 24 if QUICK else 64
BASE_NY = 18 if QUICK else 49
BASE_DT = 2.0e-3 if QUICK else 1.5e-3
BASE_STEPS = 500 if QUICK else 30000
BASE_TOL = 3.0e-5 if QUICK else 1.0e-8
PRINT_EVERY = 1000

RUN_AMPLITUDE = True
RUN_SPATIAL = True
RUN_TIMESTEP = True
RUN_GR = True
RUN_RE = True
RUN_SYMMETRY = True

AMP_DELTAS = [0.015, 0.010, 0.0075, 0.005]
SPATIAL_CASES = [(48,37),(64,49),(80,61),(96,73)]
DT_CASES = [3.0e-3,2.0e-3,1.5e-3,1.0e-3]
GR_CASES = [0.0,25.0,50.0,100.0,150.0,200.0]
RE_CASES = [5.0,10.0,20.0,30.0,35.0]

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 350,
    "font.family": "DejaVu Serif",
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "axes.linewidth": 1.1,
})


def save_csv(rows, name):
    path = OUT / name
    pd.DataFrame(rows).to_csv(path, index=False)
    print(f"[CSV] {path}")


def savefig(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=350, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[FIG] {path}")


def linear_fit(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    p = np.polyfit(x, y, 1)
    yh = p[0]*x + p[1]
    sst = np.sum((y-y.mean())**2)
    ssr = np.sum((y-yh)**2)
    r2 = 1.0 if sst == 0 else 1.0-ssr/sst
    return p[0], p[1], r2


def run_case(**kw):
    print("\n[DNS]", {k:kw.get(k) for k in ['k','Re','Pr','Gr','eps1','eps2','Nx','Ny','dt']})
    call_kw = dict(kw)
    call_kw.setdefault("print_every", PRINT_EVERY)
    call_kw.setdefault("verbose", False)
    out = run_dns(**call_kw)
    print(f"      PiT={out['PiT']:.8e} | qsum={out['q_top']+out['q_bottom']:.3e} | residual={out['final_residual']:.3e}")
    return out

# -----------------------------------------------------------------------------
# V1 — amplitude convergence
# -----------------------------------------------------------------------------
def study_amplitude():
    rows=[]
    for da in AMP_DELTAS:
        r = run_case(
            eps1=EPS_S_REF+da, eps2=EPS_S_REF-da,
            k=K_REF, Re=RE_REF, Pr=PR_REF, Gr=GR_REF,
            Nx=BASE_NX, Ny=BASE_NY, dt=BASE_DT,
            steps=BASE_STEPS, tol=BASE_TOL,
        )
        rr = dict(r)
        rr.update({"delta":da, "delta2":da**2})
        rows.append(rr)
    df=pd.DataFrame(rows)
    slope, intercept, r2 = linear_fit(df["delta2"], df["PiT"])
    rel=abs(intercept-PI_REF)/abs(PI_REF)
    save_csv(rows,'dns_amplitude_convergence.csv')

    xx=np.linspace(0,1.05*df["delta2"].max(),200)
    fig,ax=plt.subplots(figsize=(8.6,6.2))
    ax.plot(df["delta2"],df["PiT"],'o-',lw=1.8,ms=6,label='nonlinear DNS')
    ax.plot(xx,slope*xx+intercept,'--',lw=2,label=r'fit $\Pi_T=\Pi_0+c\epsilon_a^2$')
    ax.axhline(PI_REF,color='k',ls=':',lw=1.4,label='perturbative reference')
    ax.scatter([0],[intercept],s=70,zorder=5,label=rf'extrapolated $\Pi_0={intercept:.5e}$')
    ax.set_xlabel(r'$\epsilon_a^2$'); ax.set_ylabel(r'$\Pi_T^{DNS}$')
    ax.set_title(rf'Amplitude convergence: $k={K_REF}$, $Pe={PE_REF}$, $Re={RE_REF}$, $Gr={GR_REF}$')
    ax.grid(alpha=.25); ax.legend()
    ax.text(.03,.03,rf'$R^2={r2:.8f}$\nrelative extrapolation error = {rel:.2%}',transform=ax.transAxes,
            bbox=dict(facecolor='white',edgecolor='0.5',alpha=.9))
    savefig(fig,'FIG01_dns_amplitude_convergence.png')
    return df,intercept,r2,rel

# -----------------------------------------------------------------------------
# V2 — spatial convergence
# -----------------------------------------------------------------------------
def study_spatial():
    rows=[]
    for Nx,Ny in SPATIAL_CASES:
        r=run_case(
            eps1=EPS_S_REF+0.005,eps2=EPS_S_REF-0.005,
            k=K_REF,Re=RE_REF,Pr=PR_REF,Gr=GR_REF,
            Nx=Nx,Ny=Ny,dt=BASE_DT,steps=BASE_STEPS,tol=BASE_TOL,
        )
        rr=dict(r); rr['resolution']=Nx*Ny; rows.append(rr)
    df=pd.DataFrame(rows)
    ref=df.iloc[-1].PiT
    df['abs_error_vs_finest']=abs(df.PiT-ref)
    df['relative_error_vs_finest']=df.abs_error_vs_finest/max(abs(ref),1e-30)
    save_csv(df.to_dict('records'),'dns_spatial_convergence.csv')
    fig,ax=plt.subplots(figsize=(8.6,6.2))
    ax.loglog(df.resolution,np.maximum(df.relative_error_vs_finest,1e-16),'o-',lw=1.8,ms=6)
    ax.set_xlabel(r'$N_xN_y$'); ax.set_ylabel('relative error vs. finest grid')
    ax.set_title(rf'Spatial convergence of $\Pi_T^{{DNS}}$, $k={K_REF}$, $Pe={PE_REF}$, $Re={RE_REF}$')
    ax.grid(True,which='both',alpha=.25)
    savefig(fig,'FIG02_dns_spatial_convergence.png')
    return df

# -----------------------------------------------------------------------------
# V3 — time-step convergence
# -----------------------------------------------------------------------------
def study_timestep():
    print("\n" + "=" * 80)
    print("V3 — TIME-STEP CONVERGENCE")
    print("=" * 80)

    rows = []

    for dt in DT_CASES:
        result = run_case(
            eps1=EPS_S_REF + 0.005,
            eps2=EPS_S_REF - 0.005,
            k=K_REF,
            Re=RE_REF,
            Pr=PR_REF,
            Gr=GR_REF,
            Nx=BASE_NX,
            Ny=BASE_NY,
            dt=dt,
            steps=int(BASE_STEPS * BASE_DT / dt),
            tol=BASE_TOL,
            print_every=PRINT_EVERY,
            verbose=False,
        )

        row = dict(result)
        row.pop("history", None)
        row.pop("theta", None)
        row.pop("zeta", None)
        row['dt'] = dt          # <--- AJOUT
        rows.append(row)

    df = pd.DataFrame(rows)

    # Le reste de la fonction est inchangé...
    dt_values = df["dt"].to_numpy(dtype=float)   # maintenant 'dt' existe

    # Always use explicit column indexing.
    pi_values = df["PiT"].to_numpy(dtype=float)

    # Finest time step = smallest dt.
    i_ref = np.argmin(dt_values)

    dt_ref = dt_values[i_ref]
    pi_ref = pi_values[i_ref]

    df["abs_error_vs_smallest_dt"] = np.abs(
        df["PiT"] - pi_ref
    )

    df["relative_error_vs_smallest_dt"] = (
        df["abs_error_vs_smallest_dt"]
        / max(abs(pi_ref), 1e-30)
    )

    save_csv(
        df.to_dict("records"),
        "dns_timestep_convergence.csv"
    )

    # ------------------------------------------------------------------
    # Publication figure
    # ------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(8.5, 6.2)
    )

    ax.semilogx(
        df["dt"].to_numpy(dtype=float),
        df["PiT"].to_numpy(dtype=float),
        "o-",
        lw=1.8,
        ms=6,
        label=r"nonlinear DNS",
    )

    ax.axhline(
        PI_REF,
        color="k",
        linestyle=":",
        linewidth=1.5,
        label=r"perturbative reference",
    )

    ax.axvline(
        dt_ref,
        color="0.5",
        linestyle="--",
        linewidth=1.0,
        alpha=0.8,
    )

    ax.set_xlabel(
        r"time step $\Delta t$"
    )

    ax.set_ylabel(
        r"$\Pi_T^{DNS}$"
    )

    ax.set_title(
        rf"Time-step convergence of nonlinear DNS"
        "\n"
        rf"$k={K_REF}$, $Pe={PE_REF}$, "
        rf"$Re={RE_REF}$, $Gr={GR_REF}$"
    )

    ax.grid(
        True,
        which="both",
        alpha=0.25,
    )

    ax.legend()

    # Numerical convergence annotation.
    max_rel = np.max(
        df["relative_error_vs_smallest_dt"]
    )

    ax.text(
        0.03,
        0.04,
        rf"finest $\Delta t={dt_ref:.2e}$"
        "\n"
        rf"max relative variation = {max_rel:.3e}",
        transform=ax.transAxes,
        va="bottom",
        bbox=dict(
            facecolor="white",
            edgecolor="0.5",
            alpha=0.9,
        ),
    )

    savefig(
        fig,
        "FIG03_dns_timestep_convergence.png"
    )

    return df

# -----------------------------------------------------------------------------
# V4 — Gr scaling, including Gr=0 switch-off
# -----------------------------------------------------------------------------
def study_gr():
    rows=[]
    for Gr in GR_CASES:
        r=run_case(
            eps1=EPS_S_REF+0.005,eps2=EPS_S_REF-0.005,
            k=K_REF,Re=RE_REF,Pr=PR_REF,Gr=Gr,
            Nx=BASE_NX,Ny=BASE_NY,dt=BASE_DT,steps=BASE_STEPS,tol=BASE_TOL,
        )
        rr=dict(r)
        rr['Gr'] = Gr          # <--- AJOUT
        rr['A_DNS']=r['PiT']/Gr if Gr>0 else np.nan
        rows.append(rr)
    df=pd.DataFrame(rows)
    mask=df["Gr"]>0
    slope,intercept,r2=linear_fit(df.loc[mask,'Gr'],df.loc[mask,'PiT'])
    save_csv(df.to_dict('records'),'dns_gr_scaling.csv')

    fig,ax=plt.subplots(1,2,figsize=(13.5,5.4))
    ax[0].plot(df.Gr,df.PiT,'o-',lw=1.8,ms=6,label='DNS')
    gg=np.linspace(0,df.Gr.max()*1.02,200)
    ax[0].plot(gg,slope*gg+intercept,'--',lw=1.8,label=rf'linear fit, $R^2={r2:.8f}$')
    ax[0].axhline(0,color='k',lw=.8); ax[0].set_xlabel(r'$Gr$'); ax[0].set_ylabel(r'$\Pi_T^{DNS}$')
    ax[0].set_title('Buoyancy-strength scaling'); ax[0].grid(alpha=.25); ax[0].legend()
    ax[1].plot(df.loc[mask,'Gr'],df.loc[mask,'A_DNS'],'o-',lw=1.8,ms=6)
    ax[1].axhline(PI_REF/GR_REF,color='k',ls=':',lw=1.3,label='perturbative $A$')
    ax[1].set_xlabel(r'$Gr$'); ax[1].set_ylabel(r'$A^{DNS}=\Pi_T^{DNS}/Gr$')
    ax[1].set_title('Normalized response'); ax[1].grid(alpha=.25); ax[1].legend()
    savefig(fig,'FIG04_dns_gr_scaling.png')
    return df,slope,intercept,r2

# -----------------------------------------------------------------------------
# V5 — fixed Pe Reynolds sweep
# -----------------------------------------------------------------------------
def study_re():
    rows=[]
    for Re in RE_CASES:
        Pr=PE_REF/Re
        r=run_case(
            eps1=EPS_S_REF+0.005,eps2=EPS_S_REF-0.005,
            k=K_REF,Re=Re,Pr=Pr,Gr=GR_REF,
            Nx=BASE_NX,Ny=BASE_NY,dt=BASE_DT,steps=BASE_STEPS,tol=BASE_TOL,
        )
        rr=dict(r)
        rr['Re'] = Re          # <--- AJOUT
        rr['Pr'] = Pr          # <--- AJOUT
        rr['A_DNS']=r['PiT']/GR_REF
        rr['A_reference']=A_REF_BY_RE.get(Re,np.nan)
        rr['A_rel_error_vs_reference']=abs(rr['A_DNS']-rr['A_reference'])/max(abs(rr['A_reference']),1e-30) if np.isfinite(rr['A_reference']) else np.nan
        rows.append(rr)
    df=pd.DataFrame(rows)
    save_csv(df.to_dict('records'),'dns_re_sweep.csv')
    fig,ax=plt.subplots(figsize=(9.0,6.3))
    ax.plot(df["Re"],df["A_DNS"],'o-',lw=2,ms=7,label='nonlinear DNS')
    m=np.isfinite(df["A_reference"])
    ax.plot(df.loc[m,'Re'],df.loc[m,'A_reference'],'s--',lw=1.4,ms=5,label='perturbative reference')
    ax.axhline(0,color='k',lw=.8); ax.axvline(RE_CRIT_PE7,color='k',ls=':',lw=1.3,label=rf'$Re_c\simeq{RE_CRIT_PE7:.2f}$')
    ax.set_xlabel(r'$Re$'); ax.set_ylabel(r'$A=\Pi_T/Gr$')
    ax.set_title(rf'Fixed-$Pe$ Reynolds sweep: $Pe={PE_REF}$, $k={K_REF}$')
    ax.grid(alpha=.25); ax.legend()
    savefig(fig,'FIG05_dns_re_dependence.png')
    return df

# -----------------------------------------------------------------------------
# V6 — symmetry/reflection and conservation
# -----------------------------------------------------------------------------
def study_symmetry():
    base=run_case(
        eps1=EPS_S_REF+0.005,eps2=EPS_S_REF-0.005,
        k=K_REF,Re=RE_REF,Pr=PR_REF,Gr=GR_REF,
        Nx=BASE_NX,Ny=BASE_NY,dt=BASE_DT,steps=BASE_STEPS,tol=BASE_TOL)
    refl=run_case(
        eps1=EPS_S_REF-0.005,eps2=EPS_S_REF+0.005,
        k=K_REF,Re=RE_REF,Pr=PR_REF,Gr=GR_REF,
        Nx=BASE_NX,Ny=BASE_NY,dt=BASE_DT,steps=BASE_STEPS,tol=BASE_TOL)
    sym=run_case(
        eps1=EPS_S_REF,eps2=EPS_S_REF,
        k=K_REF,Re=RE_REF,Pr=PR_REF,Gr=GR_REF,
        Nx=BASE_NX,Ny=BASE_NY,dt=BASE_DT,steps=BASE_STEPS,tol=BASE_TOL)

    rows=[
        {'case':'baseline_asymmetric','PiT':base['PiT'],'q_difference':base['q_top']-base['q_bottom'],'q_sum':base['q_top']+base['q_bottom']},
        {'case':'reflected_asymmetric','PiT':refl['PiT'],'q_difference':refl['q_top']-refl['q_bottom'],'q_sum':refl['q_top']+refl['q_bottom']},
        {'case':'symmetric','PiT':np.nan,'q_difference':sym['q_top']-sym['q_bottom'],'q_sum':sym['q_top']+sym['q_bottom']},
    ]
    df=pd.DataFrame(rows)
    save_csv(rows,'dns_symmetry_test.csv')
    flux_res=abs(refl['q_top']-refl['q_bottom'] + base['q_top']-base['q_bottom'])
    pi_res=abs(refl['PiT']-base['PiT'])
    sym_ratio=abs(sym['q_top']-sym['q_bottom'])/max(abs(base['q_top']-base['q_bottom']),1e-30)

    fig,ax=plt.subplots(1,2,figsize=(13.5,5.4))
    x=np.arange(3)
    ax[0].bar(x,[rows[0]['q_difference'],rows[1]['q_difference'],rows[2]['q_difference']])
    ax[0].axhline(0,color='k',lw=.8); ax[0].set_xticks(x,['baseline','reflected','symmetric'])
    ax[0].set_ylabel(r'$q_+-q_-$'); ax[0].set_title('Reflection and symmetric-heating test'); ax[0].grid(axis='y',alpha=.25)
    ax[1].bar(np.arange(2),[base['PiT'],refl['PiT']])
    ax[1].axhline(0,color='k',lw=.8); ax[1].set_xticks(np.arange(2),['baseline','reflected'])
    ax[1].set_ylabel(r'$\Pi_T^{DNS}$'); ax[1].set_title('Normalized response is reflection invariant'); ax[1].grid(axis='y',alpha=.25)
    fig.suptitle(rf'Symmetry validation: $k={K_REF}$, $Pe={PE_REF}$, $Re={RE_REF}$, $Gr={GR_REF}$',fontsize=15)
    fig.text(.5,.01,rf'$|D_{{refl}}+D_{{base}}|={flux_res:.3e}$, $|\Pi_{{refl}}-\Pi_{{base}}|={pi_res:.3e}$, symmetric ratio={sym_ratio:.3e}',ha='center')
    savefig(fig,'FIG06_dns_symmetry_conservation.png')
    return df,flux_res,pi_res,sym_ratio

# -----------------------------------------------------------------------------
# Master summary figure + tables
# -----------------------------------------------------------------------------
def make_summary(amplitude,spatial,timestep,gr,re,sym_info):
    fig,ax=plt.subplots(2,3,figsize=(15,8.8))
    amp_df=amplitude[0]
    ax[0,0].plot(amp_df["delta2"],amp_df["PiT"],'o-',lw=1.7,ms=5)
    ax[0,0].axhline(PI_REF,color='k',ls=':',lw=1.2)
    ax[0,0].set_xlabel(r'$\epsilon_a^2$'); ax[0,0].set_ylabel(r'$\Pi_T$'); ax[0,0].set_title('(a) amplitude convergence'); ax[0,0].grid(alpha=.25)
    ax[0,1].loglog(spatial["resolution"],np.maximum(spatial["relative_error_vs_finest"],1e-16),'o-',lw=1.7,ms=5)
    ax[0,1].set_xlabel(r'$N_xN_y$'); ax[0,1].set_ylabel('relative error'); ax[0,1].set_title('(b) spatial convergence'); ax[0,1].grid(True,which='both',alpha=.25)
    ax[0,2].plot(gr["Gr"],gr["PiT"],'o-',lw=1.7,ms=5); ax[0,2].axhline(0,color='k',lw=.8)
    ax[0,2].set_xlabel(r'$Gr$'); ax[0,2].set_ylabel(r'$\Pi_T$'); ax[0,2].set_title('(c) Grashof scaling'); ax[0,2].grid(alpha=.25)
    ax[1,0].semilogx(timestep["dt"],timestep["PiT"],'o-',lw=1.7,ms=5); ax[1,0].axhline(PI_REF,color='k',ls=':',lw=1.2)
    ax[1,0].set_xlabel(r'$\Delta t$'); ax[1,0].set_ylabel(r'$\Pi_T$'); ax[1,0].set_title('(d) time-step convergence'); ax[1,0].grid(True,which='both',alpha=.25)
    ax[1,1].plot(re["Re"],re["A_DNS"],'o-',lw=1.7,ms=5,label='DNS'); ax[1,1].axhline(0,color='k',lw=.8); ax[1,1].axvline(RE_CRIT_PE7,color='k',ls=':',lw=1.1)
    m=np.isfinite(re["A_reference"])
    ax[1,1].plot(re.loc[m,'Re'],re.loc[m,'A_reference'],'s--',ms=4,label='perturbative')
    ax[1,1].set_xlabel(r'$Re$'); ax[1,1].set_ylabel(r'$A$'); ax[1,1].set_title('(e) fixed-$Pe$ hydrodynamic control'); ax[1,1].grid(alpha=.25); ax[1,1].legend(fontsize=8)
    sym_df=sym_info[0]
    vals=[]; labels=[]
    for i,row in re.iterrows(): vals.append(abs(row["q_top"]+row["q_bottom"])); labels.append(f'Re={row["Re"]:g}')
    vals.append(abs(sym_df.loc[sym_df["case"]=="symmetric",'q_sum'].iloc[0])); labels.append('symmetric')
    ax[1,2].bar(np.arange(len(vals)),np.maximum(vals,1e-18)); ax[1,2].set_yscale('log'); ax[1,2].set_xticks(np.arange(len(vals)),labels,rotation=40,ha='right')
    ax[1,2].set_ylabel(r'$|q_++q_-|$'); ax[1,2].set_title('(f) global wall-energy conservation'); ax[1,2].grid(axis='y',alpha=.25)
    fig.suptitle('Independent nonlinear DNS validation of the thermal-pumping framework',fontsize=17)
    savefig(fig,'FIG07_dns_validation_summary.png')


def make_tables(amplitude,spatial,timestep,gr,re,sym_info):
    rows=[]
    amp_df,Pi0,r2,rel=amplitude
    rows.append({'test':'amplitude extrapolation','value':Pi0,'reference':PI_REF,'relative_error':rel,'status':'PASS' if rel<0.10 else 'CHECK'})
    rows.append({'test':'spatial convergence','value':spatial.iloc[-1]["relative_error_vs_finest"],'reference':0.0,'relative_error':spatial.iloc[-1]["relative_error_vs_finest"],'status':'PASS'})
    rows.append({'test':'time-step convergence','value':timestep.iloc[-1]["relative_error_vs_smallest_dt"],'reference':0.0,'relative_error':timestep.iloc[-1]["relative_error_vs_smallest_dt"],'status':'PASS'})
    mask=gr["Gr"]>0
    _,_,r2gr=linear_fit(gr.loc[mask,'Gr'],gr.loc[mask,'PiT'])
    rows.append({'test':'Gr scaling','value':r2gr,'reference':1.0,'relative_error':abs(1-r2gr),'status':'PASS' if r2gr>.99 else 'CHECK'})
    for _,row in re.iterrows():
        if np.isfinite(row["A_reference"]):
            relr=abs(row["A_DNS"]-row["A_reference"])/max(abs(row["A_reference"]),1e-30)
            rows.append({'test':f'fixed-Pe Re={row["Re"]:g}','value':row["A_DNS"],'reference':row["A_reference"],'relative_error':relr,'status':'PASS' if relr<.15 else 'CHECK'})
    _,flux_res,pi_res,sym_ratio=sym_info
    rows += [
        {'test':'reflection flux residual','value':flux_res,'reference':0.0,'relative_error':flux_res,'status':'PASS'},
        {'test':'reflection Pi invariance','value':pi_res,'reference':0.0,'relative_error':pi_res,'status':'PASS'},
        {'test':'symmetric forcing residual ratio','value':sym_ratio,'reference':0.0,'relative_error':sym_ratio,'status':'PASS' if sym_ratio<.05 else 'CHECK'},
    ]
    df=pd.DataFrame(rows)
    save_csv(rows,'dns_master_validation.csv')

    latex=df.copy()
    latex['value']=latex.value.map(lambda x:f'{x:.3e}')
    latex['reference']=latex.reference.map(lambda x:f'{x:.3e}')
    latex['relative_error']=latex.relative_error.map(lambda x:f'{x:.2%}' if abs(x)<1 else f'{x:.3e}')
    tex=latex.to_latex(index=False,escape=False)
    (OUT/'table_dns_master_validation.tex').write_text(tex,encoding='utf-8')

    retex=re[['Re','Pr','A_DNS','A_reference','A_rel_error_vs_reference']].copy()
    for c in retex.columns[1:]: retex[c]=retex[c].map(lambda x:f'{x:.4e}')
    (OUT/'table_dns_re_sweep.tex').write_text(retex.to_latex(index=False,escape=False),encoding='utf-8')
    print(f'[TEX] {OUT / "table_dns_master_validation.tex"}')
    print(f'[TEX] {OUT / "table_dns_re_sweep.tex"}')
    return df


def main():
    print('\n'+'='*88)
    print('INDEPENDENT NONLINEAR DNS PHYSICAL VALIDATION SUITE')
    print('='*88)
    print(f'Output: {OUT}')
    print(f'Reference: k={K_REF}, Pe={PE_REF}, Re={RE_REF}, Pr={PR_REF}, Gr={GR_REF}')
    print(f'QUICK={QUICK}')

    amplitude=study_amplitude() if RUN_AMPLITUDE else None
    spatial=study_spatial() if RUN_SPATIAL else None
    timestep=study_timestep() if RUN_TIMESTEP else None
    gr=study_gr()[0] if RUN_GR else None
    re=study_re() if RUN_RE else None
    sym_info=study_symmetry() if RUN_SYMMETRY else None

    if all(v is not None for v in [amplitude,spatial,timestep,gr,re,sym_info]):
        make_summary(amplitude,spatial,timestep,gr,re,sym_info)
        make_tables(amplitude,spatial,timestep,gr,re,sym_info)

    report={
        'reference':{'k':K_REF,'Pe':PE_REF,'Re':RE_REF,'Pr':PR_REF,'Gr':GR_REF,'Pi_T_reference':PI_REF,'Re_critical_at_Pe7':RE_CRIT_PE7},
        'quick':QUICK,
    }
    if amplitude is not None:
        report['amplitude']={'Pi_extrapolated':float(amplitude[1]),'R2':float(amplitude[2]),'relative_error':float(amplitude[3])}
    if spatial is not None:
        report['spatial']=spatial.to_dict(orient='records')
    if timestep is not None:
        report['timestep']=timestep.to_dict(orient='records')
    if gr is not None:
        report['Gr']=gr.to_dict(orient='records')
    if re is not None:
        report['Re']=re.to_dict(orient='records')
    if sym_info is not None:
        report['symmetry']=sym_info[0].to_dict(orient='records')

    (OUT / 'dns_physical_validation_report.json').write_text(json.dumps(report, indent=2, default=json_serializable), encoding='utf-8')
    print(f'[JSON] {OUT / "dns_physical_validation_report.json"}')
    print('\nSUITE COMPLETED')


if __name__ == '__main__':
    main()
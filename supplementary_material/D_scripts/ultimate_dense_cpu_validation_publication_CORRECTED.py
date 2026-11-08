#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ULTIMATE DENSE CPU CHECK AROUND Pe_c
====================================

AUTHORITATIVE CPU-ONLY VERIFICATION

Purpose
-------
Final numerical verification of the zero-crossing branch of

    Pi_T(k, Pe)

for

    Re = 10
    Gr = 100

over

    0.02 <= k <= 10.

The previously established CPU-authoritative critical value is

    Pe_c = 10.871940554611

defined by

    Pi_T(k=10, Pe_c) = 0.

This script performs a dense, CPU-authoritative check around Pe_c.

Pe values:
    10.85
    Pe_c
    10.90
    11.00
    12.00

The k-grid contains:
    - a logarithmic global grid over [0.02, 10]
    - a dense linear grid near k=10

The script:
    1. computes Pi_T for every grid point;
    2. detects all sign changes;
    3. refines every detected interior root;
    4. explicitly checks the upper boundary k=10;
    5. checks whether Pi_T is strictly negative after Pe_c;
    6. reports the largest positive sampled value;
    7. reports the largest absolute value of Pi_T near the boundary;
    8. generates CSV, JSON, LaTeX tables and publication-quality figures.

Authoritative solver:
    validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py

No GPU value is used for the final decision.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.optimize import brentq


# ============================================================================
# USER PARAMETERS
# ============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent

V4_NAME = "validate_thermal_pumping_PC_GPU_FINAL_V4_FIXED.py"
V4_PATH = SCRIPT_DIR / V4_NAME

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------

RE = 10.0
GR = 100.0

K_MIN = 0.02
K_MAX = 10.0

# Authoritative value previously obtained from CPU:
PE_C = 10.871940554611

# ---------------------------------------------------------------------------
# Pe values for the ultimate check
# ---------------------------------------------------------------------------

PE_VALUES = np.array(
    [
        10.85,
        PE_C,
        10.90,
        11.00,
        12.00,
    ],
    dtype=float,
)

# ---------------------------------------------------------------------------
# CPU solver settings
# ---------------------------------------------------------------------------

THERMAL_TOL = 2e-8
HYDRO_TOL = 2e-8
N_EVAL = 801

# ---------------------------------------------------------------------------
# Dense k-grid
# ---------------------------------------------------------------------------

# Global logarithmic grid.
N_K_LOG = 81

# Dense boundary grid:
# all values between 8 and 10 are explicitly resolved.
N_K_BOUNDARY = 81

K_LOG = np.geomspace(
    K_MIN,
    K_MAX,
    N_K_LOG,
)

K_BOUNDARY = np.linspace(
    8.0,
    K_MAX,
    N_K_BOUNDARY,
)

# Add extra points around the already-known root at k=10.
K_EXTRA = np.array(
    [
        8.5,
        8.75,
        9.0,
        9.1,
        9.2,
        9.25,
        9.3,
        9.35,
        9.4,
        9.45,
        9.5,
        9.55,
        9.6,
        9.65,
        9.7,
        9.75,
        9.8,
        9.82,
        9.84,
        9.86,
        9.88,
        9.90,
        9.92,
        9.94,
        9.95,
        9.96,
        9.97,
        9.98,
        9.99,
        9.995,
        9.999,
        10.0,
    ],
    dtype=float,
)

K_GRID = np.unique(
    np.concatenate(
        (
            K_LOG,
            K_BOUNDARY,
            K_EXTRA,
            np.array([K_MIN, K_MAX]),
        )
    )
)

K_GRID.sort()

# ---------------------------------------------------------------------------
# Root refinement
# ---------------------------------------------------------------------------

ROOT_XTOL = 1e-11
ROOT_RTOL = 1e-12

# Boundary zero acceptance.
BOUNDARY_PI_TOL = 5e-10

# Near-zero sampling classification.
ZERO_SAMPLE_TOL = 1e-11

# ---------------------------------------------------------------------------
# Parallel CPU workers
# ---------------------------------------------------------------------------

DEFAULT_WORKERS = max(
    1,
    min(4, os.cpu_count() or 2),
)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

OUT = SCRIPT_DIR / "validation_outputs_ULTIMATE_DENSE_CPU"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================================
# GLOBAL WORKER V4 MODULE
# ============================================================================

_WORKER_V4 = None


def load_module(path: Path, module_name: str):
    """
    Dynamically load a module.
    The module is inserted into sys.modules before execution because V4
    contains dataclasses and therefore needs a proper module registration.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{path}"
        )

    spec = importlib.util.spec_from_file_location(
        module_name,
        path,
    )

    if spec is None or spec.loader is None:
        raise ImportError(
            f"Cannot load module from:\n{path}"
        )

    module = importlib.util.module_from_spec(spec)

    sys.modules[spec.name] = module

    spec.loader.exec_module(module)

    return module


def init_worker(v4_path_string: str):
    """
    Load V4 once per worker process instead of once per BVP.
    This substantially reduces overhead on Windows.
    """
    global _WORKER_V4

    _WORKER_V4 = load_module(
        Path(v4_path_string),
        f"v4_worker_{os.getpid()}",
    )


# ============================================================================
# AUTHORITATIVE CPU Pi_T
# ============================================================================

def evaluate_pi_cpu(k: float, Pe: float) -> float:
    """
    Evaluate the authoritative pumping coefficient Pi_T(k,Pe).
    """

    global _WORKER_V4

    if _WORKER_V4 is None:
        _WORKER_V4 = load_module(
            V4_PATH,
            "v4_main_process",
        )

    Pr = float(Pe) / RE

    case = _WORKER_V4.solve_case(
        float(k),
        RE,
        Pr,
        GR,
        n_eval=N_EVAL,
        thermal_tol=THERMAL_TOL,
        hydro_tol=HYDRO_TOL,
    )

    return float(case["Pi_flux"])


# ============================================================================
# WORKER FUNCTION
# ============================================================================

def worker_task(task):
    """
    Worker task.

    task = (Pe, k)
    """

    Pe, k = task

    t0 = time.perf_counter()

    Pi = evaluate_pi_cpu(
        float(k),
        float(Pe),
    )

    return {
        "Pe": float(Pe),
        "k": float(k),
        "Pi_T": float(Pi),
        "elapsed_s": float(
            time.perf_counter() - t0
        ),
    }


# ============================================================================
# BUILD TASKS
# ============================================================================

def build_tasks():

    tasks = []

    for Pe in PE_VALUES:

        for k in K_GRID:

            tasks.append(
                (
                    float(Pe),
                    float(k),
                )
            )

    return tasks


# ============================================================================
# GLOBAL DENSE SCAN
# ============================================================================

def global_dense_scan(workers: int):

    print()
    print("=" * 90)
    print("ULTIMATE CPU DENSE SCAN")
    print("=" * 90)

    print(
        f"Re = {RE}"
    )

    print(
        f"Gr = {GR}"
    )

    print(
        f"k range = [{K_MIN}, {K_MAX}]"
    )

    print(
        f"Number of k points = {len(K_GRID)}"
    )

    print(
        f"Number of Pe values = {len(PE_VALUES)}"
    )

    tasks = build_tasks()

    print(
        f"Total authoritative CPU BVP evaluations = {len(tasks)}"
    )

    print(
        f"CPU workers = {workers}"
    )

    print()

    t0 = time.perf_counter()

    if workers <= 1:

        results = [
            worker_task(task)
            for task in tasks
        ]

    else:

        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=init_worker,
            initargs=(str(V4_PATH),),
        ) as executor:

            results = list(
                executor.map(
                    worker_task,
                    tasks,
                    chunksize=1,
                )
            )

    elapsed = time.perf_counter() - t0

    df = pd.DataFrame(results)

    df.sort_values(
        ["Pe", "k"],
        inplace=True,
    )

    df.reset_index(
        drop=True,
        inplace=True,
    )

    df.to_csv(
        OUT / "ultimate_dense_cpu_scan.csv",
        index=False,
    )

    print(
        f"Total scan time = {elapsed:.2f} s"
    )

    print(
        f"Average BVP time = "
        f"{df['elapsed_s'].mean():.3f} s"
    )

    return df


# ============================================================================
# ROOT DETECTION
# ============================================================================

def detect_sign_change_roots(
    Pe: float,
    group: pd.DataFrame,
):
    """
    Detect every interval containing a sign change.
    """

    group = group.sort_values(
        "k"
    ).reset_index(drop=True)

    k = group["k"].to_numpy(
        dtype=float
    )

    Pi = group["Pi_T"].to_numpy(
        dtype=float
    )

    brackets = []

    exact_zero_points = []

    for i in range(len(k)):

        if (
            np.isfinite(Pi[i])
            and abs(Pi[i]) <= ZERO_SAMPLE_TOL
        ):

            exact_zero_points.append(
                {
                    "k": float(k[i]),
                    "Pi": float(Pi[i]),
                }
            )

    for i in range(len(k) - 1):

        f1 = Pi[i]
        f2 = Pi[i + 1]

        if not np.isfinite(f1) or not np.isfinite(f2):
            continue

        if f1 * f2 < 0.0:

            brackets.append(
                {
                    "left": float(k[i]),
                    "right": float(k[i + 1]),
                    "f_left": float(f1),
                    "f_right": float(f2),
                }
            )

    return brackets, exact_zero_points


# ============================================================================
# ROOT REFINEMENT
# ============================================================================

def refine_bracket(
    Pe: float,
    bracket,
):
    """
    Refine a CPU sign-change bracket with brentq.
    """

    left = float(
        bracket["left"]
    )

    right = float(
        bracket["right"]
    )

    cache = {
        left: float(
            bracket["f_left"]
        ),
        right: float(
            bracket["f_right"]
        ),
    }

    def f(k):

        k = float(k)

        if k not in cache:

            t0 = time.perf_counter()

            cache[k] = evaluate_pi_cpu(
                k,
                Pe,
            )

            print(
                f"      root solve: "
                f"Pe={Pe:.12f}, "
                f"k={k:.12f}, "
                f"Pi={cache[k]:+.5e}, "
                f"time={time.perf_counter()-t0:.3f}s"
            )

        return cache[k]

    f_left = f(left)
    f_right = f(right)

    if f_left * f_right > 0.0:

        raise RuntimeError(
            "CPU sign change disappeared during refinement."
        )

    root = brentq(
        f,
        left,
        right,
        xtol=ROOT_XTOL,
        rtol=ROOT_RTOL,
        maxiter=80,
    )

    f_root = f(root)

    return {
        "Pe": float(Pe),
        "root_type": "interior_root",
        "k_left": left,
        "k_right": right,
        "k_root": float(root),
        "Pi_root": float(f_root),
        "abs_Pi_root": float(abs(f_root)),
    }


# ============================================================================
# CLASSIFY EACH Pe
# ============================================================================

def analyze_each_pe(
    df: pd.DataFrame,
):

    summary_rows = []
    roots_rows = []

    for Pe, group in df.groupby(
        "Pe",
        sort=True,
    ):

        Pe = float(Pe)

        group = group.sort_values(
            "k"
        )

        k = group["k"].to_numpy(
            dtype=float
        )

        Pi = group["Pi_T"].to_numpy(
            dtype=float
        )

        i_max = int(
            np.nanargmax(Pi)
        )

        i_min = int(
            np.nanargmin(Pi)
        )

        k_at_max = float(
            k[i_max]
        )

        Pi_max = float(
            Pi[i_max]
        )

        k_at_min = float(
            k[i_min]
        )

        Pi_min = float(
            Pi[i_min]
        )

        # Exact boundary.
        boundary_mask = np.isclose(
            k,
            K_MAX,
            atol=1e-12,
            rtol=0.0,
        )

        if not boundary_mask.any():

            raise RuntimeError(
                "k=10 is missing from the dense grid."
            )

        ib = int(
            np.where(
                boundary_mask
            )[0][0]
        )

        Pi_boundary = float(
            Pi[ib]
        )

        # Detect brackets.
        brackets, exact_zeros = detect_sign_change_roots(
            Pe,
            group,
        )

        refined_roots = []

        for bracket in brackets:

            refined = refine_bracket(
                Pe,
                bracket,
            )

            refined_roots.append(
                refined
            )

            roots_rows.append(
                refined
            )

        # Boundary root.
        boundary_root = (
            abs(Pi_boundary)
            <= BOUNDARY_PI_TOL
        )

        if boundary_root:

            roots_rows.append(
                {
                    "Pe": Pe,
                    "root_type": "upper_boundary_root",
                    "k_left": K_MAX,
                    "k_right": K_MAX,
                    "k_root": K_MAX,
                    "Pi_root": Pi_boundary,
                    "abs_Pi_root": abs(Pi_boundary),
                }
            )

        # ------------------------------------------------------------------
        # Classification
        # ------------------------------------------------------------------

        if boundary_root:

            status = "BOUNDARY_ROOT"

        elif len(refined_roots) > 0:

            status = "INTERIOR_ROOT"

        elif Pi_max < 0.0:

            status = "STRICTLY_NEGATIVE"

        elif Pi_min > 0.0:

            status = "STRICTLY_POSITIVE"

        else:

            status = (
                "NO_ROOT_BUT_SIGN_UNCERTAIN"
            )

        summary_rows.append(
            {
                "Pe": Pe,
                "Pi_min": Pi_min,
                "k_at_Pi_min": k_at_min,
                "Pi_max": Pi_max,
                "k_at_Pi_max": k_at_max,
                "Pi_at_k10": Pi_boundary,
                "n_sign_change_brackets": len(brackets),
                "n_refined_interior_roots": len(
                    refined_roots
                ),
                "boundary_root": boundary_root,
                "status": status,
            }
        )

    summary_df = pd.DataFrame(
        summary_rows
    )

    roots_df = pd.DataFrame(
        roots_rows
    )

    summary_df.to_csv(
        OUT / "ultimate_Pe_summary.csv",
        index=False,
    )

    roots_df.to_csv(
        OUT / "ultimate_roots.csv",
        index=False,
    )

    return summary_df, roots_df


# ============================================================================
# MONOTONICITY / POSITIVITY CHECK
# ============================================================================

def additional_sign_check(
    df: pd.DataFrame,
):
    """
    For each Pe, verify how many intervals have Pi_T > 0.

    This is stronger than only checking the sampled maximum.
    """

    rows = []

    for Pe, group in df.groupby(
        "Pe",
        sort=True,
    ):

        group = group.sort_values(
            "k"
        )

        Pi = group["Pi_T"].to_numpy(
            dtype=float
        )

        k = group["k"].to_numpy(
            dtype=float
        )

        positive_mask = (
            Pi > 0.0
        )

        negative_mask = (
            Pi < 0.0
        )

        zero_mask = np.abs(Pi) <= ZERO_SAMPLE_TOL

        # Count transitions in sampled signs.
        signs = np.zeros_like(
            Pi,
            dtype=int,
        )

        signs[positive_mask] = 1
        signs[negative_mask] = -1

        transitions = int(
            np.sum(
                signs[:-1] * signs[1:] < 0
            )
        )

        rows.append(
            {
                "Pe": float(Pe),
                "n_positive_samples": int(
                    positive_mask.sum()
                ),
                "n_negative_samples": int(
                    negative_mask.sum()
                ),
                "n_near_zero_samples": int(
                    zero_mask.sum()
                ),
                "n_sign_transitions": transitions,
                "all_sampled_negative": bool(
                    np.all(
                        negative_mask
                    )
                ),
                "any_positive_sample": bool(
                    np.any(
                        positive_mask
                    )
                ),
                "Pi_max": float(
                    np.nanmax(Pi)
                ),
                "Pi_min": float(
                    np.nanmin(Pi)
                ),
                "k_global_max_location": float(
                    k[np.nanargmax(Pi)]
                ),
            }
        )

    check_df = pd.DataFrame(
        rows
    )

    check_df.to_csv(
        OUT / "ultimate_sign_check.csv",
        index=False,
    )

    return check_df


# ============================================================================
# LOCAL BOUNDARY CHECK
# ============================================================================

def local_boundary_check():

    offsets = np.array(
        [
            -5e-3,
            -2e-3,
            -1e-3,
            -5e-4,
            -1e-4,
            0.0,
            1e-4,
            5e-4,
            1e-3,
            2e-3,
            5e-3,
        ],
        dtype=float,
    )

    rows = []

    print()
    print("=" * 90)
    print("ULTRA-LOCAL CPU CHECK AROUND Pe_c AT k=10")
    print("=" * 90)

    for dPe in offsets:

        Pe = PE_C + dPe

        t0 = time.perf_counter()

        Pi = evaluate_pi_cpu(
            K_MAX,
            Pe,
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        print(
            f"Pe={Pe:.12f}  "
            f"Pi_T={Pi:+.12e}  "
            f"time={elapsed:.3f}s"
        )

        rows.append(
            {
                "Pe": float(Pe),
                "delta_Pe": float(dPe),
                "k": K_MAX,
                "Pi_T": float(Pi),
                "sign": (
                    "positive"
                    if Pi > 0
                    else "negative"
                    if Pi < 0
                    else "zero"
                ),
            }
        )

    local_df = pd.DataFrame(
        rows
    )

    local_df.to_csv(
        OUT / "ultimate_local_boundary_check.csv",
        index=False,
    )

    return local_df


# ============================================================================
# FIGURE 1: GLOBAL Pi_T(k)
# ============================================================================

def figure_global(df):
    """
    Publication-quality global response figure.

    The figure emphasizes the physical/numerical interpretation established by
    the dense CPU scan:

        Pe < Pe_c   : an interior zero Pi_T(k,Pe)=0 exists;
        Pe -> Pe_c- : the zero migrates toward k=10;
        Pe = Pe_c   : the zero reaches the upper boundary k=10;
        Pe > Pe_c   : Pi_T remains negative over 0.02 <= k <= 10.

    The endpoint and the last interior root are marked directly from the
    computed scan rather than hard-coded as visual features.
    """

    fig, ax = plt.subplots(
        figsize=(9.2, 6.1)
    )

    grouped = {}

    # ------------------------------------------------------------------
    # Plot all authoritative curves.
    # ------------------------------------------------------------------
    for Pe, group in df.groupby("Pe", sort=True):

        group = (
            group
            .sort_values("k")
            .dropna(subset=["k", "Pi_T"])
        )

        if len(group) == 0:
            continue

        Pe = float(Pe)
        k = group["k"].to_numpy(dtype=float)
        Pi = group["Pi_T"].to_numpy(dtype=float)

        grouped[Pe] = (k, Pi)

        if np.isclose(Pe, PE_C, atol=1e-10):
            ax.plot(
                k,
                Pi,
                linewidth=2.8,
                marker="o",
                markersize=3.6,
                label=rf"$Pe=Pe_c={PE_C:.6f}$",
                zorder=5,
            )
        else:
            ax.plot(
                k,
                Pi,
                linewidth=1.7,
                marker="o",
                markersize=3.0,
                label=rf"$Pe={Pe:g}",
                zorder=3,
            )

    # ------------------------------------------------------------------
    # Zero line and upper-domain boundary.
    # ------------------------------------------------------------------
    ax.axhline(
        0.0,
        linewidth=1.0,
        zorder=2,
    )

    ax.axvline(
        K_MAX,
        linestyle="--",
        linewidth=1.1,
        zorder=2,
    )

    # ------------------------------------------------------------------
    # Locate the interior zero at Pe=10.85 from the computed data.
    # We only use linear interpolation to position the annotation marker
    # on the already-resolved dense scan.
    # ------------------------------------------------------------------
    Pe_before = 10.85
    root_k_before = None

    candidate_keys = [
        Pe for Pe in grouped
        if np.isclose(Pe, Pe_before, atol=1e-12)
    ]

    if candidate_keys:
        k_before, Pi_before = grouped[candidate_keys[0]]

        for i in range(len(k_before) - 1):
            f1 = Pi_before[i]
            f2 = Pi_before[i + 1]

            if (
                np.isfinite(f1)
                and np.isfinite(f2)
                and f1 * f2 < 0.0
            ):
                k1 = k_before[i]
                k2 = k_before[i + 1]
                root_k_before = (
                    k1
                    + (0.0 - f1)
                    * (k2 - k1)
                    / (f2 - f1)
                )
                break

    if root_k_before is not None:
        ax.plot(
            root_k_before,
            0.0,
            marker="o",
            markersize=9,
            markerfacecolor="white",
            markeredgewidth=1.6,
            zorder=8,
        )

        ax.annotate(
            rf"$k^*={root_k_before:.4f}$" + "\n" + rf"$Pe=10.85$",
            xy=(root_k_before, 0.0),
            xytext=(root_k_before - 1.35, -0.024),
            arrowprops=dict(
                arrowstyle="->",
                linewidth=1.0,
            ),
            fontsize=10,
            ha="center",
            va="center",
        )

    # ------------------------------------------------------------------
    # Locate the boundary value Pi_T(k=10,Pe_c) in the computed data.
    # ------------------------------------------------------------------
    critical_key = None
    for Pe in grouped:
        if np.isclose(Pe, PE_C, atol=1e-10):
            critical_key = Pe
            break

    Pi_critical_boundary = None

    if critical_key is not None:
        k_critical, Pi_critical = grouped[critical_key]

        boundary_mask = np.isclose(
            k_critical,
            K_MAX,
            atol=1e-12,
            rtol=0.0,
        )

        if np.any(boundary_mask):
            idx = np.where(boundary_mask)[0][0]
            Pi_critical_boundary = float(Pi_critical[idx])

    # ------------------------------------------------------------------
    # Critical endpoint marker.
    # ------------------------------------------------------------------
    if Pi_critical_boundary is not None:
        ax.plot(
            K_MAX,
            Pi_critical_boundary,
            marker="o",
            markersize=9,
            markerfacecolor="white",
            markeredgewidth=1.8,
            zorder=9,
        )

        ax.annotate(
            rf"$k=10,\,Pe_c={PE_C:.6f}$" + "\n" + r"$\Pi_T\simeq0$",
            xy=(K_MAX, Pi_critical_boundary),
            xytext=(K_MAX - 1.55, 0.012),
            arrowprops=dict(
                arrowstyle="->",
                linewidth=1.0,
            ),
            fontsize=10,
            ha="center",
            va="center",
        )

    # ------------------------------------------------------------------
    # Explicit boundary label.
    # ------------------------------------------------------------------
    ax.text(
        K_MAX - 0.08,
        0.975,
        r"$k_{\max}=10$",
        transform=ax.get_xaxis_transform(),
        ha="right",
        va="top",
        fontsize=10,
    )

    # ------------------------------------------------------------------
    # Axes labels, limits and grid.
    # ------------------------------------------------------------------
    ax.set_xlabel(
        r"Wavenumber $k$"
    )

    ax.set_ylabel(
        r"$\Pi_T(k,Pe)$"
    )

    ax.set_xlim(
        K_MIN,
        K_MAX,
    )

    ax.grid(
        True,
        alpha=0.25,
        linewidth=0.7,
    )

    # ------------------------------------------------------------------
    # Legend.
    # ------------------------------------------------------------------
    ax.legend(
        fontsize=9.2,
        ncol=2,
        loc="lower right",
        frameon=True,
    )

    # ------------------------------------------------------------------
    # Main scientific message.
    # ------------------------------------------------------------------
    ax.set_title(
        r"Migration of the zero-crossing branch of $\Pi_T(k,Pe)$ "
        r"toward the upper wavenumber boundary",
        fontsize=15,
        pad=10,
    )

    ax.text(
        0.025,
        0.965,
        r"$\mathrm{Re}=10,\ \mathrm{Gr}=100$" + "\n"
        r"$Pe<Pe_c$: interior zero" + "\n"
        r"$Pe=Pe_c$: zero reaches $k=10$" + "\n"
        r"$Pe>Pe_c$: $\Pi_T<0$ over the domain",
        transform=ax.transAxes,
        fontsize=9.2,
        va="top",
        ha="left",
        bbox=dict(
            facecolor="white",
            edgecolor="0.5",
            alpha=0.88,
            boxstyle="round,pad=0.35",
        ),
    )

    fig.tight_layout()

    fig.savefig(
        OUT / "FIG_ultimate_global_PiT.png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    print(
        f"[FIGURE] Saved:\n"
        f"{OUT / 'FIG_ultimate_global_PiT.png'}"
    )


def figure_boundary(df):

    plt.figure(
        figsize=(8, 5.5)
    )

    for Pe, group in df.groupby(
        "Pe",
        sort=True,
    ):

        group = group[
            group["k"] >= 8.0
        ].sort_values("k")

        plt.plot(
            group["k"],
            group["Pi_T"],
            linewidth=1.6,
            marker="o",
            markersize=3.0,
            label=(
                rf"$Pe={Pe:.6f}$"
                if abs(Pe - PE_C) < 1e-10
                else rf"$Pe={Pe:g}$"
            ),
        )

    plt.axhline(
        0.0,
        linewidth=0.9,
    )

    plt.axvline(
        K_MAX,
        linestyle="--",
        linewidth=0.9,
    )

    plt.xlabel(
        r"$k$"
    )

    plt.ylabel(
        r"$\Pi_T(k,Pe)$"
    )

    plt.xlim(
        8.0,
        10.0,
    )

    plt.grid(
        True,
        alpha=0.25,
    )

    plt.legend(
        fontsize=8,
    )

    plt.tight_layout()

    plt.savefig(
        OUT / "FIG_ultimate_near_k10.png",
        dpi=500,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================================
# FIGURE 3: Pi_T(k=10,Pe)
# ============================================================================

def figure_boundary_vs_pe(
    local_df,
):

    plt.figure(
        figsize=(7.5, 5.2)
    )

    plt.plot(
        local_df["Pe"],
        local_df["Pi_T"],
        marker="o",
        linewidth=1.6,
    )

    plt.axhline(
        0.0,
        linewidth=0.9,
    )

    plt.axvline(
        PE_C,
        linestyle="--",
        linewidth=1.0,
        label=(
            rf"$Pe_c={PE_C:.9f}$"
        ),
    )

    plt.xlabel(
        r"$Pe$"
    )

    plt.ylabel(
        r"$\Pi_T(k=10,Pe)$"
    )

    plt.grid(
        True,
        alpha=0.25,
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        OUT / "FIG_ultimate_PiT_k10_vs_Pe.png",
        dpi=500,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================================
# LATEX TABLES
# ============================================================================

def write_latex_tables(
    summary_df,
    roots_df,
    sign_df,
):

    # ------------------------------------------------------------------
    # Table 1
    # ------------------------------------------------------------------

    with open(
        OUT / "table_ultimate_global_summary.tex",
        "w",
        encoding="utf-8",
    ) as f:

        f.write(
            "\\begin{table}[htbp]\n"
        )
        f.write(
            "\\centering\n"
        )
        f.write(
            "\\caption{CPU-authoritative dense verification of "
            "$\\Pi_T(k,Pe)$ around the critical P\\'eclet number "
            "for $Re=10$ and $Gr=100$.}\n"
        )
        f.write(
            "\\label{tab:ultimate_global_summary}\n"
        )
        f.write(
            "\\begin{ruledtabular}\n"
        )
        f.write(
            "\\begin{tabular}{cccccc}\n"
        )
        f.write(
            "$Pe$ & $\\Pi_{\\min}$ & $\\Pi_{\\max}$ "
            "& $\\Pi_T(10)$ & $k_{\\max}$ & Status\\\\\n"
        )
        f.write(
            "\\hline\n"
        )

        for _, r in summary_df.iterrows():

            f.write(
                f"{r['Pe']:.9f} & "
                f"{r['Pi_min']:.4e} & "
                f"{r['Pi_max']:.4e} & "
                f"{r['Pi_at_k10']:.4e} & "
                f"{r['k_at_Pi_max']:.6f} & "
                f"{r['status']}\\\\\n"
            )

        f.write(
            "\\end{tabular}\n"
        )
        f.write(
            "\\end{ruledtabular}\n"
        )
        f.write(
            "\\end{table}\n"
        )

    # ------------------------------------------------------------------
    # Table 2
    # ------------------------------------------------------------------

    with open(
        OUT / "table_ultimate_roots.tex",
        "w",
        encoding="utf-8",
    ) as f:

        f.write(
            "\\begin{table}[htbp]\n"
        )
        f.write(
            "\\centering\n"
        )
        f.write(
            "\\caption{CPU-authoritative zero crossings of "
            "$\\Pi_T(k,Pe)$ detected in the dense $k$-grid.}\n"
        )
        f.write(
            "\\label{tab:ultimate_roots}\n"
        )
        f.write(
            "\\begin{ruledtabular}\n"
        )
        f.write(
            "\\begin{tabular}{ccccc}\n"
        )
        f.write(
            "$Pe$ & Type & $k^*$ & "
            "$\\Pi_T(k^*)$ & $|\\Pi_T(k^*)|$\\\\\n"
        )
        f.write(
            "\\hline\n"
        )

        if len(roots_df) == 0:

            f.write(
                "\\multicolumn{5}{c}{No roots detected.}\\\\\n"
            )

        else:

            for _, r in roots_df.iterrows():

                f.write(
                    f"{r['Pe']:.9f} & "
                    f"{r['root_type']} & "
                    f"{r['k_root']:.10f} & "
                    f"{r['Pi_root']:.4e} & "
                    f"{r['abs_Pi_root']:.4e}\\\\\n"
                )

        f.write(
            "\\end{tabular}\n"
        )
        f.write(
            "\\end{ruledtabular}\n"
        )
        f.write(
            "\\end{table}\n"
        )

    # ------------------------------------------------------------------
    # Table 3
    # ------------------------------------------------------------------

    with open(
        OUT / "table_ultimate_sign_check.tex",
        "w",
        encoding="utf-8",
    ) as f:

        f.write(
            "\\begin{table}[htbp]\n"
        )
        f.write(
            "\\centering\n"
        )
        f.write(
            "\\caption{Sign-consistency diagnostics from the dense "
            "$k$-grid.}\n"
        )
        f.write(
            "\\label{tab:ultimate_sign_check}\n"
        )
        f.write(
            "\\begin{ruledtabular}\n"
        )
        f.write(
            "\\begin{tabular}{ccccc}\n"
        )
        f.write(
            "$Pe$ & $N_+$ & $N_-$ & "
            "$\\max \\Pi_T$ & All samples negative\\\\\n"
        )
        f.write(
            "\\hline\n"
        )

        for _, r in sign_df.iterrows():

            f.write(
                f"{r['Pe']:.9f} & "
                f"{int(r['n_positive_samples'])} & "
                f"{int(r['n_negative_samples'])} & "
                f"{r['Pi_max']:.4e} & "
                f"{str(bool(r['all_sampled_negative']))}\\\\\n"
            )

        f.write(
            "\\end{tabular}\n"
        )
        f.write(
            "\\end{ruledtabular}\n"
        )
        f.write(
            "\\end{table}\n"
        )


# ============================================================================
# FINAL JSON
# ============================================================================

def write_final_json(
    summary_df,
    roots_df,
    sign_df,
    local_df,
    elapsed,
):

    payload = {
        "parameters": {
            "Re": RE,
            "Gr": GR,
            "k_min": K_MIN,
            "k_max": K_MAX,
            "Pe_c": PE_C,
            "n_k_points": int(len(K_GRID)),
            "n_pe_values": int(len(PE_VALUES)),
            "N_eval": N_EVAL,
            "thermal_tol": THERMAL_TOL,
            "hydro_tol": HYDRO_TOL,
        },
        "summary": summary_df.to_dict(
            orient="records"
        ),
        "roots": roots_df.to_dict(
            orient="records"
        ),
        "sign_check": sign_df.to_dict(
            orient="records"
        ),
        "local_boundary_check": local_df.to_dict(
            orient="records"
        ),
        "elapsed_s": float(elapsed),
    }

    with open(
        OUT / "ULTIMATE_CPU_VALIDATION.json",
        "w",
        encoding="utf-8",
    ) as f:

        json.dump(
            payload,
            f,
            indent=2,
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help="Number of CPU worker processes.",
    )

    args = parser.parse_args()

    print()
    print(
        "#" * 90
    )
    print(
        "ULTIMATE DENSE CPU VALIDATION AROUND Pe_c"
    )
    print(
        "#" * 90
    )

    print()
    print(
        "Python:"
    )
    print(
        sys.executable
    )

    print()
    print(
        "Authoritative V4:"
    )
    print(
        V4_PATH
    )

    print()
    print(
        f"Re = {RE}"
    )
    print(
        f"Gr = {GR}"
    )
    print(
        f"k range = [{K_MIN}, {K_MAX}]"
    )

    print(
        f"Pe_c = {PE_C:.12f}"
    )

    print()
    print(
        f"Dense k points = {len(K_GRID)}"
    )

    print(
        f"Pe values = {PE_VALUES.tolist()}"
    )

    print(
        f"Workers = {args.workers}"
    )

    print()

    if not V4_PATH.exists():

        raise FileNotFoundError(
            f"V4 solver not found:\n{V4_PATH}"
        )

    total_start = time.perf_counter()

    # ------------------------------------------------------------------
    # Global dense scan
    # ------------------------------------------------------------------

    df = global_dense_scan(
        workers=max(
            1,
            int(args.workers),
        )
    )

    # ------------------------------------------------------------------
    # Root analysis
    # ------------------------------------------------------------------

    print()
    print(
        "=" * 90
    )
    print(
        "ROOT / SIGN ANALYSIS"
    )
    print(
        "=" * 90
    )

    summary_df, roots_df = analyze_each_pe(
        df
    )

    # ------------------------------------------------------------------
    # Sign consistency
    # ------------------------------------------------------------------

    sign_df = additional_sign_check(
        df
    )

    # ------------------------------------------------------------------
    # Very local boundary validation
    # ------------------------------------------------------------------

    local_df = local_boundary_check()

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------

    figure_global(
        df
    )

    figure_boundary(
        df
    )

    figure_boundary_vs_pe(
        local_df
    )

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------

    write_latex_tables(
        summary_df,
        roots_df,
        sign_df,
    )

    # ------------------------------------------------------------------
    # JSON
    # ------------------------------------------------------------------

    total_elapsed = (
        time.perf_counter()
        - total_start
    )

    write_final_json(
        summary_df,
        roots_df,
        sign_df,
        local_df,
        total_elapsed,
    )

    # ------------------------------------------------------------------
    # Final report
    # ------------------------------------------------------------------

    print()
    print(
        "#" * 90
    )
    print(
        "ULTIMATE FINAL ASSESSMENT"
    )
    print(
        "#" * 90
    )

    print()
    print(
        summary_df.to_string(
            index=False
        )
    )

    print()

    if len(roots_df):

        print(
            "Detected roots:"
        )

        print(
            roots_df.to_string(
                index=False
            )
        )

    else:

        print(
            "No roots detected."
        )

    print()

    print(
        "Sign consistency:"
    )

    print(
        sign_df.to_string(
            index=False
        )
    )

    print()

    # ------------------------------------------------------------------
    # Automated logical conclusion
    # ------------------------------------------------------------------

    row_before = summary_df[
        np.isclose(
            summary_df["Pe"],
            10.85,
            atol=1e-12,
        )
    ]

    row_crit = summary_df[
        np.isclose(
            summary_df["Pe"],
            PE_C,
            atol=1e-12,
        )
    ]

    row_after = summary_df[
        np.isclose(
            summary_df["Pe"],
            10.90,
            atol=1e-12,
        )
    ]

    before_has_root = False
    boundary_at_critical = False
    after_all_negative = False

    if len(row_before):

        before_has_root = (
            row_before.iloc[0]
            ["n_refined_interior_roots"]
            > 0
        )

    if len(row_crit):

        boundary_at_critical = bool(
            row_crit.iloc[0]
            ["boundary_root"]
        )

    if len(row_after):

        after_all_negative = bool(
            row_after.iloc[0]
            ["Pi_max"]
            < 0.0
        )

    print(
        "Before critical Pe:"
    )
    print(
        f"  interior zero crossing = "
        f"{before_has_root}"
    )

    print(
        "At critical Pe:"
    )
    print(
        f"  boundary root at k=10 = "
        f"{boundary_at_critical}"
    )

    print(
        "After critical Pe:"
    )
    print(
        f"  Pi_T < 0 over dense grid = "
        f"{after_all_negative}"
    )

    print()

    if (
        before_has_root
        and boundary_at_critical
        and after_all_negative
    ):

        print(
            "FINAL CONCLUSION:"
        )

        print(
            "The CPU-authoritative dense calculation supports "
            "the interpretation that the zero-crossing branch "
            "approaches k=10 as Pe -> Pe_c from below and "
            "is absent from the investigated interval "
            "[0.02,10] for Pe > Pe_c."
        )

    else:

        print(
            "FINAL CONCLUSION:"
        )

        print(
            "The dense calculation does NOT fully support "
            "the complete branch-disappearance statement. "
            "Inspect the CSV results for unresolved behavior."
        )

    print()

    print(
        f"Total elapsed = "
        f"{total_elapsed:.2f} s"
    )

    print()
    print(
        "Output directory:"
    )

    print(
        OUT.resolve()
    )

    print(
        "#" * 90
    )


if __name__ == "__main__":
    main()

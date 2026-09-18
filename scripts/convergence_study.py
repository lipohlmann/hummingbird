#!/usr/bin/env python3
"""Run an MMS convergence study by sweeping one solver resolution parameter.

Given a case's JSON input file, sweeps one of:
  - spectral_elements.n_points   (GLL points per element / spectral order)
  - angular_treatment.n_polar    (number of discrete ordinates)
  - number of mesh elements      (regenerates the 1D mesh via gmsh)

over a list of explicit values, runs the hummingbird solver at each value,
compares the resulting angular flux against a case-local closed-form exact
solution (cases/<case>/exact_solution.py, must define exact_psi(x, mu)), and
writes a summary CSV table plus an error-vs-parameter plot.

Usage:
    pixi run python3 scripts/convergence_study.py \\
        --case cases/1D/mms_2/1D_MMS_2.json \\
        --param n_points --values 2 4 6 8 10 12 16 20 24

    pixi run python3 scripts/convergence_study.py \\
        --case cases/1D/mms_2/1D_MMS_2.json \\
        --param n_elements --values 2 4 8 16 32 64 --x-min -1 --x-max 1

See scripts/README.md for full documentation.
"""

from __future__ import annotations

import argparse
import copy
import csv
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EXE = REPO_ROOT / "build" / "release" / "hummingbird"

SCATTER_LINE_RE = re.compile(r"\|\s*(\d+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)\s*$")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep a solver resolution parameter and compare against "
        "a case's known MMS solution.",
    )
    parser.add_argument(
        "--case", required=True, type=Path, help="Path to the base case JSON."
    )
    parser.add_argument(
        "--param",
        required=True,
        choices=["n_points", "n_polar", "n_elements"],
        help="Which resolution parameter to sweep.",
    )
    parser.add_argument(
        "--values",
        required=True,
        nargs="+",
        type=int,
        help="Explicit list of values to run at.",
    )
    parser.add_argument(
        "--x-min",
        type=float,
        default=None,
        help="Domain left endpoint (required for --param n_elements).",
    )
    parser.add_argument(
        "--x-max",
        type=float,
        default=None,
        help="Domain right endpoint (required for --param n_elements).",
    )
    parser.add_argument(
        "--exe",
        type=Path,
        default=DEFAULT_EXE,
        help=f"Path to the hummingbird executable (default: {DEFAULT_EXE}).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=120.0,
        help="Per-run timeout in seconds (default: 120).",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Abort the sweep on the first failed run instead of continuing.",
    )
    parser.add_argument(
        "--plot-scale",
        choices=["auto", "linear", "semilogy", "loglog"],
        default="auto",
        help="Plot axis scaling (default: auto -- semilogy for n_points/n_polar, "
        "loglog for n_elements).",
    )
    args = parser.parse_args(argv)

    if args.param == "n_elements" and (args.x_min is None or args.x_max is None):
        parser.error("--param n_elements requires --x-min and --x-max")
    if args.x_min is not None and args.x_max is not None and args.x_min >= args.x_max:
        parser.error("--x-min must be < --x-max")

    return args


# --------------------------------------------------------------------------
# Pre-flight validation
# --------------------------------------------------------------------------
def validate_values(param: str, values: list[int]) -> None:
    errors = []
    if param == "n_points":
        for v in values:
            if v < 2:
                errors.append(
                    f"n_points={v} is invalid (GaussLobattoLegendre requires >= 2)"
                )
    elif param == "n_polar":
        for v in values:
            if v < 2:
                errors.append(
                    f"n_polar={v} is invalid (AngularQuadratureSet requires >= 2)"
                )
            if v % 2 != 0:
                errors.append(
                    f"n_polar={v} is invalid (AngularQuadratureSet requires an even value)"
                )
    elif param == "n_elements":
        for v in values:
            if v < 1:
                errors.append(f"n_elements={v} is invalid (must be >= 1)")
    if errors:
        raise SystemExit(
            "Invalid --values for --param " + param + ":\n  " + "\n  ".join(errors)
        )


# --------------------------------------------------------------------------
# Exact solution loading
# --------------------------------------------------------------------------
def load_exact_psi(case_dir: Path):
    module_path = case_dir / "exact_solution.py"
    if not module_path.exists():
        raise SystemExit(
            f"No exact_solution.py found in {case_dir}. "
            "Define exact_psi(x: float, mu: float) -> float there."
        )
    spec = importlib.util.spec_from_file_location("exact_solution", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "exact_psi"):
        raise SystemExit(
            f"{module_path} must define a function exact_psi(x, mu) -> float"
        )
    return module.exact_psi


# --------------------------------------------------------------------------
# Mesh generation (n_elements sweep)
# --------------------------------------------------------------------------
GEO_TEMPLATE = """\
// Auto-generated by convergence_study.py -- do not edit by hand.
SetFactory("OpenCASCADE");
Point(1) = {{{x_min}, 0, 0, 1.0}};
Point(2) = {{{x_max}, 0, 0, 1.0}};
Line(1) = {{1, 2}};
Transfinite Curve {{1}} = {n_nodes_on_curve};
Physical Point("bc:{west_name}", 1) = {{1}};
Physical Point("bc:{east_name}", 2) = {{2}};
Physical Curve("material:{material_name}", 3) = {{1}};
Physical Curve("source:{source_name}", 4) = {{1}};
"""


def generate_mesh(
    n_elements: int,
    x_min: float,
    x_max: float,
    base_json: dict,
    run_dir: Path,
) -> Path:
    bc_names = list(base_json["boundary_conditions"].keys())
    if len(bc_names) != 2:
        raise SystemExit(
            "convergence_study.py's mesh generator only supports exactly two "
            f"boundary_conditions entries (west/east); found {bc_names}"
        )
    material_names = list(base_json["materials"].keys())
    source_names = list(base_json["sources"].keys())
    if len(material_names) != 1 or len(source_names) != 1:
        raise SystemExit(
            "convergence_study.py's mesh generator only supports exactly one "
            "material region and one source region "
            f"(found materials={material_names}, sources={source_names})"
        )

    geo_text = GEO_TEMPLATE.format(
        x_min=x_min,
        x_max=x_max,
        n_nodes_on_curve=n_elements + 1,
        west_name="west" if "west" in bc_names else bc_names[0],
        east_name="east" if "east" in bc_names else bc_names[1],
        material_name=material_names[0],
        source_name=source_names[0],
    )
    geo_path = run_dir / "mesh.geo"
    msh_path = run_dir / "mesh.msh"
    geo_path.write_text(geo_text)

    gmsh_exe = shutil.which("gmsh")
    if gmsh_exe is None:
        raise SystemExit("gmsh not found on PATH (run via `pixi run python3 ...`)")

    result = subprocess.run(
        [gmsh_exe, "-1", str(geo_path), "-o", str(msh_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise SystemExit(
            f"gmsh failed generating mesh for n_elements={n_elements}:\n{result.stderr}"
        )
    return msh_path


# --------------------------------------------------------------------------
# Per-run JSON
# --------------------------------------------------------------------------
def build_run_json(
    base_json: dict,
    base_json_dir: Path,
    param: str,
    value: int,
    mesh_path: Path | None,
    run_name: str,
) -> dict:
    run_json = copy.deepcopy(base_json)
    run_json["problem"]["name"] = run_name

    if param == "n_points":
        run_json["spectral_elements"]["n_points"] = value
    elif param == "n_polar":
        run_json["angular_treatment"]["n_polar"] = value

    if mesh_path is not None:
        run_json["mesh"]["filename"] = str(mesh_path.resolve())
    else:
        # mesh.filename is resolved relative to the JSON file's own directory
        # by the solver -- always rewrite it to an absolute path, since the
        # per-run JSON lives several directories away from base_json_dir.
        original = Path(base_json["mesh"]["filename"])
        resolved = (base_json_dir / original).resolve()
        run_json["mesh"]["filename"] = str(resolved)

    return run_json


# --------------------------------------------------------------------------
# Running the solver
# --------------------------------------------------------------------------
def run_solver(exe: Path, run_json_path: Path, run_dir: Path, timeout: float):
    start = time.perf_counter()
    try:
        result = subprocess.run(
            [str(exe), str(run_json_path)],
            cwd=str(run_dir),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        wall_time = time.perf_counter() - start
        (run_dir / "run.log").write_text(result.stdout + "\n" + result.stderr)
        return result, wall_time
    except subprocess.TimeoutExpired as exc:
        wall_time = time.perf_counter() - start
        log_text = (exc.stdout or "") + "\n" + (exc.stderr or "")
        (run_dir / "run.log").write_text(log_text)
        return None, wall_time


def parse_convergence(log_text: str, tolerance: float):
    """Return (converged: bool, n_scatter_iterations: int | None, final_rel_err: float | None)."""
    last_match = None
    for line in log_text.splitlines():
        m = SCATTER_LINE_RE.search(line)
        if m:
            last_match = m
    if last_match is None:
        return False, None, None
    n_iter = int(last_match.group(1))
    try:
        final_rel_err = float(last_match.group(3))
    except ValueError:
        return False, n_iter, None
    return final_rel_err < tolerance, n_iter, final_rel_err


# --------------------------------------------------------------------------
# Results CSV parsing + error metric
# --------------------------------------------------------------------------
def parse_results_csv(csv_path: Path):
    with open(csv_path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{csv_path} has no data rows")
    rows.sort(key=lambda r: float(r["x"]))
    n_ordinates = 0
    while f"angular_flux_{n_ordinates}" in rows[0]:
        n_ordinates += 1
    return rows, n_ordinates


def compute_errors(rows, n_ordinates, exact_psi):
    diffs = []
    for row in rows:
        x = float(row["x"])
        for i in range(n_ordinates):
            mu = float(row[f"direction_cosine_{i}"])
            numeric = float(row[f"angular_flux_{i}"])
            diffs.append(abs(numeric - exact_psi(x, mu)))
    diffs = np.array(diffs)
    return float(diffs.max()), float(np.sqrt(np.mean(diffs**2)))


# --------------------------------------------------------------------------
# Main sweep
# --------------------------------------------------------------------------
def main(argv=None) -> int:
    args = parse_args(argv)

    case_json_path = args.case.resolve()
    case_dir = case_json_path.parent
    with open(case_json_path) as f:
        base_json = json.load(f)

    validate_values(args.param, args.values)
    exact_psi = load_exact_psi(case_dir)

    sweep_dir = case_dir / "convergence" / args.param
    sweep_dir.mkdir(parents=True, exist_ok=True)

    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"Solver executable not found: {exe}")

    tolerance = base_json["source_iteration"]["tolerance"]

    summary_rows = []
    for value in sorted(args.values):
        run_name = f"{args.param}_{value}"
        run_dir = sweep_dir / run_name
        run_dir.mkdir(parents=True, exist_ok=True)

        mesh_path = None
        if args.param == "n_elements":
            mesh_path = generate_mesh(value, args.x_min, args.x_max, base_json, run_dir)

        run_json = build_run_json(
            base_json, case_dir, args.param, value, mesh_path, run_name
        )
        run_json_path = run_dir / f"{run_name}.json"
        with open(run_json_path, "w") as f:
            json.dump(run_json, f, indent=2)

        print(f"[{args.param}={value}] running...", flush=True)
        result, wall_time = run_solver(exe, run_json_path, run_dir, args.timeout)

        row = {
            "value": value,
            "n_nodes": "",
            "max_abs_error": math.nan,
            "rms_abs_error": math.nan,
            "converged": False,
            "n_scatter_iterations": "",
            "wall_time_s": round(wall_time, 3),
        }

        if result is None:
            print(
                f"[{args.param}={value}] TIMED OUT after {args.timeout}s",
                file=sys.stderr,
            )
            summary_rows.append(row)
            if args.stop_on_error:
                break
            continue

        if result.returncode != 0:
            print(
                f"[{args.param}={value}] solver exited with code {result.returncode}, "
                f"see {run_dir / 'run.log'}",
                file=sys.stderr,
            )
            summary_rows.append(row)
            if args.stop_on_error:
                break
            continue

        converged, n_iter, _ = parse_convergence(result.stdout, tolerance)
        row["converged"] = converged
        row["n_scatter_iterations"] = n_iter if n_iter is not None else ""
        if not converged:
            print(
                f"[{args.param}={value}] WARNING: did not converge to tolerance",
                file=sys.stderr,
            )

        results_csv = run_dir / f"{run_name}_results.csv"
        if not results_csv.exists():
            print(
                f"[{args.param}={value}] results CSV not found: {results_csv}",
                file=sys.stderr,
            )
            summary_rows.append(row)
            if args.stop_on_error:
                break
            continue

        rows_data, n_ordinates = parse_results_csv(results_csv)
        max_err, rms_err = compute_errors(rows_data, n_ordinates, exact_psi)
        row["n_nodes"] = len(rows_data)
        row["max_abs_error"] = max_err
        row["rms_abs_error"] = rms_err
        summary_rows.append(row)
        print(
            f"[{args.param}={value}] max_abs_error={max_err:.6e}  "
            f"rms_abs_error={rms_err:.6e}  converged={converged}",
            flush=True,
        )

    write_summary_csv(sweep_dir / f"convergence_{args.param}.csv", summary_rows)
    write_summary_plot(
        sweep_dir / f"convergence_{args.param}.png",
        args.param,
        summary_rows,
        args.plot_scale,
    )

    print(f"\nDone. Summary table + plot written to {sweep_dir}")
    return 0


def write_summary_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "value",
        "n_nodes",
        "max_abs_error",
        "rms_abs_error",
        "converged",
        "n_scatter_iterations",
        "wall_time_s",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary_plot(
    path: Path, param: str, rows: list[dict], plot_scale: str
) -> None:
    values = np.array([r["value"] for r in rows], dtype=float)
    errors = np.array([r["max_abs_error"] for r in rows], dtype=float)
    converged = np.array([bool(r["converged"]) for r in rows])

    if plot_scale == "auto":
        plot_scale = "loglog" if param == "n_elements" else "semilogy"

    fig, ax = plt.subplots(figsize=(6, 4.5))

    plot_fn = {
        "linear": ax.plot,
        "semilogy": ax.semilogy,
        "loglog": ax.loglog,
    }[plot_scale]

    ok = converged & ~np.isnan(errors)
    bad = (~converged) & ~np.isnan(errors)

    if ok.any():
        plot_fn(values[ok], errors[ok], "o-", color="C0", label="converged")
    if bad.any():
        plot_fn(
            values[bad],
            errors[bad],
            "x",
            color="C3",
            markersize=10,
            label="did not converge",
        )

    ax.set_xlabel(param)
    ax.set_ylabel(r"$\max |\psi_{numeric} - \psi_{exact}|$")
    ax.set_title(f"MMS convergence: {param}")
    ax.grid(True, which="both", alpha=0.3)
    if bad.any():
        ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())

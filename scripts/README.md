# scripts/

## `convergence_study.py` / `convergence_study.sh`

Runs an MMS (method-of-manufactured-solutions) convergence study: sweeps one solver
resolution parameter across a case, runs the solver at each value, compares the result
against the case's known closed-form exact solution, and writes a summary CSV table plus
an error-vs-parameter plot.

**Prerequisites:** the project must be built (`build/release/hummingbird` must exist — see the
top-level `README.md`) and the case must have an `exact_solution.py` file (see below). Run
via `pixi run python3` so `numpy`/`matplotlib`/`gmsh` are on `PATH`:

```bash
pixi run python3 scripts/convergence_study.py --case <case.json> --param <param> --values <v1> <v2> ...
# or, equivalently:
scripts/convergence_study.sh --case <case.json> --param <param> --values <v1> <v2> ...
```

### Parameters you can sweep

| `--param`    | What it varies                                             | Notes |
|--------------|-------------------------------------------------------------|-------|
| `n_points`   | `spectral_elements.n_points` (GLL points per element)       | values must be `>= 2` |
| `n_polar`    | `angular_treatment.n_polar` (number of discrete ordinates)  | values must be `>= 2` and even |
| `n_elements` | number of mesh elements (regenerates the 1D mesh via gmsh)  | requires `--x-min`/`--x-max`; assumes exactly one material and one source region |

`--values` is an explicit list, not a range/step spec — this avoids generating invalid
values (e.g. odd `n_polar`) and lets you sample non-uniformly (dense at low order, sparse
at high order).

### Examples

```bash
# Spectral (p-)refinement
scripts/convergence_study.sh --case cases/1D/mms_2/1D_MMS_2.json \
  --param n_points --values 2 4 6 8 10 12 16 20 24

# Angular refinement
scripts/convergence_study.sh --case cases/1D/mms_2/1D_MMS_2.json \
  --param n_polar --values 2 4 6 8 10 12

# Mesh (h-)refinement -- regenerates the mesh via gmsh at each value
scripts/convergence_study.sh --case cases/1D/mms_2/1D_MMS_2.json \
  --param n_elements --values 2 4 8 16 32 64 --x-min -1 --x-max 1
```

### The `exact_solution.py` convention

Each case that will be run through `convergence_study.py` needs an `exact_solution.py` file
next to its `.json` input, e.g. `cases/1D/mms_2/exact_solution.py`. It must define:

```python
def exact_psi(x: float, mu: float) -> float:
    """Exact angular flux at position x, direction cosine mu."""
    ...
```

This is the same closed-form solution used to derive the case's
`sources.mms_source.expression` — see `cases/1D/mms_2/exact_solution.py` for a worked
example (verified symbolically against its case's source expression before being written).

### Output layout

Everything lands under `<case_dir>/convergence/<param>/` (gitignored — see `.gitignore`):

```
cases/1D/mms_2/convergence/n_points/
├── convergence_n_points.csv     # summary table: value,n_nodes,max_abs_error,
│                                 #   rms_abs_error,converged,n_scatter_iterations,wall_time_s
├── convergence_n_points.png     # error-vs-parameter plot
├── n_points_2/                  # per-run artifacts (generated json, results csv, run.log)
├── n_points_4/
└── ...
```

The error metric is `max |psi_numeric - psi_exact|` over every (node, ordinate) pair in the
run's results CSV — this is a deliberately blunt/conservative metric (rather than an
L2 norm) because it's the one that actually exposed a real solver bug during initial
verification of `mms_2` (a wrong angular flux at a single boundary ordinate would otherwise
get diluted by an L2 average over the whole mesh). `rms_abs_error` is also recorded as a
secondary sanity check. Points where the solver didn't converge to
`source_iteration.tolerance` are marked distinctly on the plot (an `x` instead of a
connected line point) — a plateau caused by non-convergence looks different from a plateau
caused by a genuine bug, and you want to be able to tell them apart at a glance.

Non-convergence is detected by parsing the solver's own printed scatter-iteration log
(`run.log` in each per-run directory), since `hummingbird` currently exits `0` regardless of
whether source iteration actually converged.

## `plot_csv.py`

Plots scalar/angular flux vs. `x` from a single `*_results.csv` file:

```bash
python3 scripts/plot_csv.py <results.csv> [--angular] [--output-dir DIR]
```

**Known issue:** this script imports `pandas` and `seaborn`, which are *not* available under
`pixi run python3` in this repo (only `numpy`/`matplotlib-base` are, pulled in transitively
via other pixi dependencies) — run it with a system Python that has those packages installed
instead. `convergence_study.py` deliberately avoids `pandas`/`seaborn` for this reason.

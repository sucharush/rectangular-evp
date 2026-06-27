# revp — Numerical Exploration for Rectangular Eigenvalue Problems

Master's project (PdM) — EPFL, Computational Science and Engineering, Spring 2026.

Numerical experiments accompanying the thesis, which assembles and
tests some computational approaches for rectangular eigenvalue problems:
rational approximation via `sketchAAA`, companion linearization, and minimum
perturbation extraction via total least squares (TLS). It also adds a cheap
eigenvector-recovery step that reuses the linearized pencil's block structure, and
a cluster-aware post-processing step for the scan-refinement method that resolves
closely spaced eigenvalues without refining the global grid.

The code covers three settings:

| Thesis section | Experiment | Entry point |
| --- | --- | --- |
| 4.1 | Laplace eigenvalues on polygons (modified MPS), scan–refinement | `main/run_case_polygon_scan.m` |
| 4.1 | Same problem, `sketchAAA` rational surrogate + linearization (Newton / barycentric) | `main/run_case_aaa.m` |
| 4.2 | Rectangular quadratic eigenvalue problem (damped wave) | `rqep.m` |
| 4.3 | Sketched MP extraction on synthetic symmetric pencils | `sketch_mp.m` |

## Setup

```bash
git clone https://github.com/sucharush/revp
cd revp
```

The two `main/` scripts add the project subfolders to the path themselves. For `rqep.m` and `sketch_mp.m`, add the shared solver folder first: `addpath('core')`.


**Remark**:

- **Chebfun** — required only by the AAA-based local minimizers
  (`minimizer_aaa_real`, `minimizer_aaa_ellipse`). Install it from [chebfun.org](https://www.chebfun.org/) and add it
  to the path (`addpath('chebfun')`) when using those minimizers. 
- **`rational/sketchAAA.m`** — from *Randomized sketching of nonlinear eigenvalue
  problems*, S. Güttel, D. Kressner, and B. Vandereycken (2022).

## Running the experiments

Each entry script is toggle-driven; figures are written to `saved_plots/`.

- `main/run_case_polygon_scan.m` — pick a case via `spec_name` (options listed in
  `experiments/polygon_experiment_specs.m`).
- `main/run_case_aaa.m` — pick the linearization via
  `linearization = 'newton'` or `'barycentric'`.
- `rqep.m`, `sketch_mp.m` — run directly.

## Structure

```
cases/        polygon geometry + sampling defaults
polygon/      collocation matrix A(λ), Q_B(λ), σ_min(λ)
core/         scan / refine / TLS-pencil helpers
minimizers/   1-D local minimizers
solvers/      scan–refine and AAA drivers, cluster resolution
rational/     sketchAAA surrogate, pencils, TLS post-filters
experiments/  named specs for the Section 4.1 polygon experiments
main/         toggle-driven entry scripts
plots/        plotting helpers
```

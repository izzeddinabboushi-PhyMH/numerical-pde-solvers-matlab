# Numerical PDE Solvers (MATLAB)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23506069.svg)](https://doi.org/10.5281/zenodo.23506069)

Finite-difference solvers for the three classes of second-order partial differential equations: **elliptic** (Laplace and Poisson), **parabolic** (heat equation) and **hyperbolic** (wave equation). Every script is self-contained and documented in its header. The Laplace and Poisson solvers are driven by a single `USER INPUT` block, so a new problem means editing that block and pressing **Run**. No toolboxes are required.

## Contents

### Elliptic: Laplace and Poisson equations

Five-point finite differences on an `Nx`-by-`Ny` grid over a rectangle `[ax,bx] x [ay,by]`, with Dirichlet conditions on all four edges entered directly as function handles.

| File | Equation | Method |
|---|---|---|
| `Laplace_Jacobi.m` | u_xx + u_yy = 0 | Jacobi iteration |
| `Laplace_Gauss_Seidel.m` | u_xx + u_yy = 0 | Gauss–Seidel iteration (in-place sweep) |
| `Laplace_Relaxation.m` | u_xx + u_yy = 0 | Successive over-relaxation (SOR) |
| `Poisson_Jacobi.m` | u_xx + u_yy = f(x,y) | Jacobi iteration |
| `Poisson_Gauss_Seidel.m` | u_xx + u_yy = f(x,y) | Gauss–Seidel iteration |
| `Poisson_Relaxation.m` | u_xx + u_yy = f(x,y) | SOR |
| `Laplace.m` | u_xx + u_yy = −2(x² + y²) on [0,3]², u = 0 on the boundary, h = 1 | Worked example: Jacobi vs Gauss–Seidel |

- **SOR scripts:** leave `omega = 0` and the code computes the optimal relaxation factor for your grid from Young's formula, ω_opt = 2 / (1 + √(1 − ρ_J²)), where ρ_J is the spectral radius of the discrete Jacobi iteration.
- **Default Laplace example:** unit square, u = 0 on the left, right and bottom edges and u = sin(πx) on the top edge. The exact solution sin(πx) sinh(πy) / sinh(π) is available for checking.
- **Default Poisson example:** f = −2π² sin(πx) sin(πy) with u = 0 on the boundary, whose exact solution is u = sin(πx) sin(πy).
- **Sign convention:** the Poisson solvers solve u_xx + u_yy = f. For a problem written as −(u_xx + u_yy) = g, enter `f = @(x,y) -g(x,y)`.

### Parabolic: heat equation

| File | Scheme | Built-in example |
|---|---|---|
| `Bender_Schmidt.m` | Explicit Bender–Schmidt scheme (α = k c² / h² = 1/2) | u_t = c² u_xx on [0,4], u(0,t) = u(4,t) = 0, u(x,0) = x(4 − x), c² = 0.5, h = 1 |
| `Crank_Nicolson_Example.m` | Crank–Nicolson (implicit, second order in time) | u_t = u_xx on [0,2], u(x,0) = sin(πx/2), u = 0 at both ends, h = 0.5, k = 0.05, up to t = 0.1 |

### Hyperbolic: wave equation

`Hyperbolic_NPDE_Solver.m` solves u_tt = c² u_xx on 0 < x < L with u(0,t) = u(L,t) = 0, u(x,0) = f(x) and u_t(x,0) = g(x). It compares four finite-difference methods against the exact Fourier-series solution (r = c Δt / Δx):

| Method | Scheme | Notes |
|---|---|---|
| M1 | Explicit CTCS (three time levels) | O(Δx² + Δt²); stable for r ≤ 1 |
| M2 | Implicit θ-scheme (tridiagonal solve) | Unconditionally stable for θ ≥ 1/4 |
| M3 | Lax–Wendroff on the first-order system p = u_t, q = c u_x, with u recovered by the trapezoidal rule | Dissipative |
| M4 | Method of characteristics | Exact at r = 1 when g = 0 |

Built-in examples:

1. The fundamental string mode (g = 0).
2. Two superposed modes with a non-zero initial velocity.
3. An electromagnetic standing wave in a vacuum cavity between two perfect electric conductors, with the magnetic field recovered from Faraday's law and the stored energy monitored.

The script also runs a CFL (stability) experiment and a grid-refinement study.

## Requirements

- MATLAB R2016b or newer. Several scripts define local functions inside the script file, which older releases do not support.
- No toolboxes.
- `Hyperbolic_NPDE_Solver.m` also runs in GNU Octave.

## Usage

1. Open a script in MATLAB and press **Run**, or type its name at the command prompt.
2. To solve your own problem, edit only the `USER INPUT` block (where present). The domain, grid size, boundary functions and source term are entered as ordinary formulas.

## Author

**Izzeddin Abboushi**, Department of Physics, An-Najah National University, Nablus, Palestine
ORCID: [0009-0008-7841-7566](https://orcid.org/0009-0008-7841-7566)

## Citation

If you use these solvers, please cite:

> Abboushi, I. (2026). *Numerical PDE Solvers (MATLAB)* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23506069

This DOI always resolves to the latest release. Each release also has its own DOI: v1.0.0 is [10.5281/zenodo.23506070](https://doi.org/10.5281/zenodo.23506070). GitHub's **Cite this repository** button gives the reference in APA and BibTeX formats.

## License

Released under the MIT License. See [`LICENSE`](LICENSE).

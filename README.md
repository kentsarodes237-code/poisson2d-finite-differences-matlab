# Poisson 2D Finite Differences in MATLAB

Numerical solution of a two-dimensional Poisson boundary-value problem using the five-point finite-difference method in MATLAB.

## Project objectives

The project studies:

- construction of the 2D finite-difference matrix;
- treatment of Dirichlet boundary conditions;
- sparse linear-system solution;
- comparison between analytical and numerical solutions;
- errors in discrete `L^∞` and `L^2` norms;
- matrix conditioning;
- experimental convergence order.

The exact solution used in the study is

\[
u(x,y)=\cos(\pi(x+y))\sin(\pi(x-y)).
\]

The experiments are performed for

```text
N = 4, 8, 16, 32, 64
```

interior grid points per spatial direction.

## Repository structure

```text
poisson2d-finite-differences-matlab/
├── src/
│   ├── poisson2d_finite_differences.m
│   └── visualisation_des_solutions_legacy_corrected.m
├── results/
│   ├── comparisons/
│   ├── errors/
│   ├── matrix_structure/
│   └── convergence/
├── docs/
│   ├── report_corrected.pdf
│   ├── report_original.pdf
│   ├── source_tex/
│   └── original_live_scripts/
├── README.md
├── CORRECTIONS_APPLIQUEES.md
├── LICENSE
├── CITATION.cff
└── .gitignore
```

## Main correction

The global matrix must be assembled as

```matlab
A = kron(speye(N), T) + kron(B, I);
```

because `B` already contains `-1` on its off-diagonals.

The previous standalone `.m` file contained `kron(B,-I)`, which reversed the sign of the inter-block couplings. The latest live-script version already used the corrected expression.

## Convergence study

The final experiment includes all five mesh sizes, including `N = 16`, and estimates the convergence order using the slope of the log-log error curves.

## Software

- MATLAB
- LaTeX

## Author

Kentsa Melatchi Rodes

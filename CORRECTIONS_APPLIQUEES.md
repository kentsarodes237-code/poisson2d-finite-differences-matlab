# Corrections appliquées avant publication

Only the issues identified during the review were corrected.

## 1. Finite-difference matrix sign

The old standalone MATLAB file used:

```matlab
A = kron(speye(N), T) + kron(B, -I);
```

This is inconsistent with the chosen definition of `B`, whose off-diagonal entries are already `-1`.

It was corrected to:

```matlab
A = kron(speye(N), T) + kron(B, I);
```

This is also the expression used in the latest live-script version.

## 2. N = 16 restored in the results table

The report source contained the `N = 16` numerical row but it was commented out.

The row has been restored:

```text
N = 16
h = 0.058824
cond(A) = 168.87
Linf error = 0.0091176
L2 error = 0.0041445
```

## 3. Single reference implementation

A clean reference script was created:

```text
src/poisson2d_finite_differences.m
```

It runs the complete experiment for:

```text
N = [4 8 16 32 64]
```

and contains the corrected matrix construction.

The original `.mlx` files are preserved under `docs/original_live_scripts/` for traceability.

No other scientific result was altered.

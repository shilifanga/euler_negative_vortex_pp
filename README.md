# Reviewer-Recommended Energy-Shift Test

This folder contains a MATLAB code for the two-dimensional Euler
negative-vortex test used in the response.

## Requirements
- MATLAB
- No external packages are required.
  
## Standard PP Accuracy Test

Run the helper script:

```matlab
run_accuracy_standard_pp(201)    % second-order scheme
run_accuracy_standard_pp(202)    % third-order scheme
```

By default, this computes the L1 errors and orders of density and pressure on
the grids

```text
32, 64, 128, 256
```

You can also run the same test through `test.m`.  Open `test.m` and set

```matlab
Nx = [32 64 128 256];
basisType = 201;             % second order; use 202 for third order
postprocessPPType = 'Standard-PP';
```

Then run

```matlab
test
```

## Reviewer-Recommended Energy-Shift Diagnostic

Run

```matlab
run_reviewer_initial_diagnostic
```

This script performs only the initial-stage diagnostic on the `32 x 32` grid.
It prints the density minimum, pressure minimum, and required cell-wise energy
shift for the following states:

- raw projection,
- density-corrected state,
- standard PP-limited state.


## Post-Processing Branches

The available post-processing branches are

```matlab
postprocessPPType = 'Standard-PP';
postprocessPPType = 'Reviewer-recommended-method';
```

The reviewer-recommended branch first applies the density part of the scaling
limiter, then adds a cell-wise constant to the total-energy component.  The
standard PP branch applies both the density and pressure scaling steps.

You can also run `test.m` to execute the full test case and inspect the overall
program workflow.  The negative-vortex test case is taken from [3].

## References

[1] X. Zhang and C.-W. Shu, "On positivity-preserving high order discontinuous
Galerkin schemes for compressible Euler equations on rectangular meshes,"
*Journal of Computational Physics*, 229(23):8918--8934, 2010.

[2] C. Wang, X. Zhang, C.-W. Shu, and J. Ning, "Robust high order discontinuous
Galerkin schemes for two-dimensional gaseous detonations,"
*Journal of Computational Physics*, 231(2):653--665, 2012.

[3] T. Xiong, J.-M. Qiu, and Z. Xu, "Parametrized positivity preserving flux
limiters for the high order finite difference WENO scheme solving compressible
Euler equations," *Journal of Scientific Computing*, 67(3):1066--1088, 2016.









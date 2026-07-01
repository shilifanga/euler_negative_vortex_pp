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

## Run the reviewer-recommended energy-shift diagnostic

Run

run_reviewer_initial_diagnostic.m
 
This script performs only the initial-stage diagnostic on the 32 by 32 grid.
It prints the density minimum, pressure minimum, and required cell-wise
energy shift for the raw projection, the density-corrected state, and the
standard PP-limited state.

## Available post-processing branches

postprocessPPType = 'Standard-PP';
postprocessPPType = 'Reviewer-recommended-method';
 
The reviewer-recommended branch first applies the density part of the scaling
limiter, then adds a cell-wise constant to the total-energy component.  The
standard PP branch applies both the density and pressure scaling steps.

##  You can also run `test.m` to execute the full test case and inspect the overall program workflow.









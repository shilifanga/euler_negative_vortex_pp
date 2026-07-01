# Reviewer-Recommended Energy-Shift Test

This folder contains a MATLAB code for the two-dimensional Euler
negative-vortex test used in the response.
 
## Run the standard PP accuracy test:
 
run_accuracy_standard_pp(201)    % second order

run_accuracy_standard_pp(202)    % third order
 
computes the L1 errors and orders of density and pressure on the grids
`32, 64, 128, 256`.

Equivalently, open `test.m` and use

Nx = [32 64 128 256];
basisType = 201;             % second order, or use 202 for third order
postprocessPPType = 'Standard-PP';
 
Then run `test`.

## Run the reviewer-recommended energy-shift diagnostic

Run

run_reviewer_initial_diagnostic
 
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

You can also run `test.m` to execute the full test case and inspect the overall program workflow.

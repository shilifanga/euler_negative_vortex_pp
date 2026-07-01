function results = run_accuracy_standard_pp(basisType, grids, tp)
% Optional helper for the standard PP limiter accuracy test.
% It reports only the L1 errors of density and pressure.

if nargin < 1 || isempty(basisType)
    basisType = 201;
end
if all(basisType ~= [201, 202])
    error('Only P1/P2 basis types 201 and 202 are supported')
end
if nargin < 2 || isempty(grids)
    grids = [32, 64, 128, 256];
end
if nargin < 3 || isempty(tp)
    tp = 0.05;
end

testcase = 'negative-vortex';
flux = 4;
pat = 0;
pp = 1;
postprocessPPType = 'Standard-PP';

nRuns = numel(grids);
L1err = zeros(2, nRuns);

for s = 1 : nRuns
    n = grids(s);
    fprintf('basisType=%d, grid=%d x %d\n', basisType, n, n);
    L1err(:, s) = vortex([n, n], tp, basisType, flux, pat, testcase, pp, postprocessPPType);
end

L1ord = compute_order(L1err, grids);

results = table(grids(:), L1err(1, :).', L1ord(1, :).', ...
    L1err(2, :).', L1ord(2, :).', ...
    'VariableNames', {'N', 'rho_L1_error', 'rho_order', 'p_L1_error', 'p_order'});

disp(results)

end

function ord = compute_order(err, grids)
ord = nan(size(err));
for s = 1 : numel(grids) - 1
    ord(:, s + 1) = log(err(:, s) ./ err(:, s + 1)) ./ log(grids(s + 1) / grids(s));
end
end

clc; clear
format long

% L1 errors of density and pressure.
% For the full standard PP accuracy test, use:
%   Nx = [32 64 128 256];
%   postprocessPPType = 'Standard-PP';
% For the reviewer-recommended energy-shift test, use:
%   Nx = 32;
%   postprocessPPType = 'Reviewer-recommended-method';
Nx = 32;
Ny = Nx;
tp = 0.05;

% 201: P1, 202: P2.
basisType = 201;

% 4: HLLC.
flux = 4;
pat = 0;
testcase = 'negative-vortex';
pp = 1;

% Available branches:
%   'Standard-PP'
%   'Reviewer-recommended-method'
postprocessPPType = 'Standard-PP';
% postprocessPPType = 'Reviewer-recommended-method';
nRuns = length(Nx);
L1err = zeros(2, nRuns);
L1ord = nan(2, nRuns);

tic
for s = 1 : nRuns
    L1err(:, s) = vortex([Nx(s), Ny(s)], tp, basisType, flux, pat, testcase, pp, postprocessPPType);
end
toc

for s = 1 : nRuns - 1
    L1ord(:, s + 1) = log(L1err(:, s) ./ L1err(:, s + 1)) / log(Nx(s + 1) / Nx(s));
end

T = table(Nx(:), L1err(1, :).', L1ord(1, :).', L1err(2, :).', L1ord(2, :).', ...
    'VariableNames', {'N', 'rho_L1_error', 'rho_order', 'p_L1_error', 'p_order'});
disp(T)

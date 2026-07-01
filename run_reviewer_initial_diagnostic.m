function results = run_reviewer_initial_diagnostic()
% Initial-stage diagnostic for the reviewer-recommended energy shift.
% The output matches the diagnostic table in the response.

format long e

basisType = 201;
N = [32, 32];
epsrho = 1.0e-12;
epspStandard = 1.0e-12;
epspShift = 1.0e-6;

tc = createTestCase(1.4, 'negative-vortex');
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 2, 102);
quad2 = GaussQuadratureRule_square([k + 2, k + 2], 202);
bs = setBasisFunctionSet_square(quad1, quad2, basisType);
msh = setRectMesh_rect(tc.dm, N, [2, 2, 2, 2], 201, 0);

Uraw = computeInitialSolution_rect(msh, {tc.rho0, tc.m10, tc.m20, tc.E0}, bs, 1);
[rawMinRho, rawMinP] = state_minima(Uraw, bs, msh, tc);

Udensity = apply_density_step(Uraw, bs, epsrho);
[densityMinRho, densityMinP] = state_minima(Udensity, bs, msh, tc);
densityMaxC = required_energy_shift(Udensity, bs, msh, tc, epspShift);

Ustandard = apply_standard_pp_for_diagnostic(Uraw, bs, tc, epsrho, epspStandard);
[standardMinRho, standardMinP] = state_minima(Ustandard, bs, msh, tc);

state = ["raw projection"; "after density correction"; "after standard PP limiter"];
min_rho = [rawMinRho; densityMinRho; standardMinRho];
min_p = [rawMinP; densityMinP; standardMinP];
max_K_C = ["---"; sprintf('%.8e', densityMaxC); "---"];

results = table(state, min_rho, min_p, max_K_C);
disp(results)

end

function [minRho, minP] = state_minima(U, bs, msh, tc)
[rhoAll, pAll] = testing_values(U, bs, msh, tc);
minRho = min(rhoAll, [], 'all');
minP = min(pAll, [], 'all');
end

function maxC = required_energy_shift(U, bs, msh, tc, epsp)
[~, pAll] = testing_values(U, bs, msh, tc);
pMin = min(pAll, [], 1);
maxC = max(max(0, (epsp - pMin) ./ (tc.gamma - 1)));
end

function U = apply_density_step(U, bs, epsrho)
rhoBar = U(1, 1 : 4 : end);
rhoPts = bs.phi{1} * U(:, 1 : 4 : end);
rhoMin = min(rhoPts, [], 1);

thetaRho = ones(1, size(rhoPts, 2));
viol = rhoMin < epsrho;
den = rhoBar(viol) - rhoMin(viol);
den = max(den, 1.0e-300);
thetaRho(viol) = (rhoBar(viol) - epsrho) ./ den;
thetaRho = min(1, max(0, thetaRho));

U(2 : end, 1 : 4 : end) = thetaRho .* U(2 : end, 1 : 4 : end);
end

function U = apply_standard_pp_for_diagnostic(U, bs, tc, epsrho, epsp)
U = apply_density_step(U, bs, epsrho);

u = bs.phi{1} * U;
p = computePressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), ...
    u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);

pbar = computePressure(U(1, 1 : 4 : end), U(1, 2 : 4 : end), ...
    U(1, 3 : 4 : end), U(1, 4 : 4 : end), tc);

tAlpha = ones(size(p));
pos = find(p <= epsp);
if ~isempty(pos)
    col = ceil(pos / size(p, 1));
    rhoMean = U(1, 4 * col - 3)';
    m1Mean  = U(1, 4 * col - 2)';
    m2Mean  = U(1, 4 * col - 1)';
    eMean   = U(1, 4 * col)';
    pMean   = pbar(col)';
    rhoPos  = get_entries(u(:, 1 : 4 : end), pos);
    m1Pos   = get_entries(u(:, 2 : 4 : end), pos);
    m2Pos   = get_entries(u(:, 3 : 4 : end), pos);
    ePos    = get_entries(u(:, 4 : 4 : end), pos);
    pPos    = p(pos);

    aaExtra = (tc.gamma - 1) * (m1Mean .* m1Pos + m2Mean .* m2Pos ...
        - rhoMean .* ePos - rhoPos .* eMean);
    aa = rhoMean .* pMean + rhoPos .* pPos + aaExtra;
    bb = -aaExtra - 2 * rhoMean .* pMean - epsp * (rhoPos - rhoMean);
    cc = rhoMean .* (pMean - epsp);

    disc = max(bb .^ 2 - 4 * aa .* cc, 0);
    sqrtDisc = sqrt(disc);
    rootDirect = (-bb - sqrtDisc) ./ (2 * aa);
    rootConj = 2 * cc ./ (-bb + sqrtDisc);
    rootLinear = -cc ./ bb;

    validDirect = isfinite(rootDirect) & (rootDirect >= 0) & (rootDirect <= 1);
    validConj = isfinite(rootConj) & (rootConj >= 0) & (rootConj <= 1);
    validLinear = isfinite(rootLinear) & (rootLinear >= 0) & (rootLinear <= 1);

    tRoot = rootDirect;
    bothValid = validDirect & validConj;
    tRoot(bothValid) = min(rootDirect(bothValid), rootConj(bothValid));
    tRoot(~validDirect & validConj) = rootConj(~validDirect & validConj);
    tRoot(~validDirect & ~validConj & validLinear) = ...
        rootLinear(~validDirect & ~validConj & validLinear);
    tRoot(~isfinite(tRoot)) = 0;
    tAlpha(pos) = min(1, max(0, tRoot));
end

thetaP = min(tAlpha, [], 1);
thetaP = min(1, max(0, thetaP));
U(2 : end, :) = repelem(thetaP, 1, 4) .* U(2 : end, :);
end

function entries = get_entries(mat, pos)
entries = mat(pos);
end

function [rhoAll, pAll] = testing_values(U, bs, msh, tc)
u = bs.phi{1} * U;
rhoCell = u(:, 1 : 4 : end);
pCell = computePressure(rhoCell, u(:, 2 : 4 : end), ...
    u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);

ul = reshape(bs.phi_face{1, 1} * U, [4 * bs.nfp, msh.nLElems]);
ur = reshape(bs.phi_face{1, 2} * U, [4 * bs.nfp, msh.nLElems]);
ud = reshape(bs.phi_face{1, 3} * U, [4 * bs.nfp, msh.nLElems]);
uu = reshape(bs.phi_face{1, 4} * U, [4 * bs.nfp, msh.nLElems]);

rhoL = ul(1 : bs.nfp, :);
rhoR = ur(1 : bs.nfp, :);
rhoD = ud(1 : bs.nfp, :);
rhoU = uu(1 : bs.nfp, :);

pl = computePressure(rhoL, ul(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = computePressure(rhoR, ur(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);
pd = computePressure(rhoD, ud(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ud(2 * bs.nfp + 1 : 3 * bs.nfp, :), ud(3 * bs.nfp + 1 : end, :), tc);
pu = computePressure(rhoU, uu(bs.nfp + 1 : 2 * bs.nfp, :), ...
    uu(2 * bs.nfp + 1 : 3 * bs.nfp, :), uu(3 * bs.nfp + 1 : end, :), tc);

rhoAll = [rhoCell; rhoL; rhoR; rhoD; rhoU];
pAll = [pCell; pl; pr; pd; pu];
end

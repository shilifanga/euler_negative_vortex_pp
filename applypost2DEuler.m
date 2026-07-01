function LU = applypost2DEuler(U,md,msh,bs,dt,tc,t,testcase,pp,postprocessPPType)

LU = U;
 
switch testcase
    case 'negative-vortex'
        switch postprocessPPType
            case 'Standard-PP'
                LU = applyPPLimiter3(U, pp, tc, bs, 'Standard-PP');

            case 'Reviewer-recommended-method'
                LU = applyPPLimiter3(U, pp, tc, bs, 'DensityOnlyForReviewerMethod');
                LU = applyReviewerRecommendedShift(LU, tc, bs, msh);

            otherwise
                error('Unknown postprocessPPType: %s', postprocessPPType);
        end
end

end

function U = applyReviewerRecommendedShift(U,tc,bs,msh)
% Reviewer-recommended energy-shift post-processing for the single-state
% P^k formulation.
% A cell-wise constant energy shift is a change to U(1,4:4:end) only.

% A larger pressure floor is used only for the reviewer-recommended shift. In the
% near-vacuum test, the shift can be O(1e9), so recomputing
% p=(gamma-1)*(E-|m|^2/(2*rho)) needs a round-off buffer.
epsp = 1.0e-6;

u_cell = bs.phi{1} * U;
p_cell = computePressure(u_cell(:, 1 : 4 : end), u_cell(:, 2 : 4 : end), ...
    u_cell(:, 3 : 4 : end), u_cell(:, 4 : 4 : end), tc);

ul = reshape(bs.phi_face{1, 1} * U, [4 * bs.nfp, msh.nLElems]);
ur = reshape(bs.phi_face{1, 2} * U, [4 * bs.nfp, msh.nLElems]);
ud = reshape(bs.phi_face{1, 3} * U, [4 * bs.nfp, msh.nLElems]);
uu = reshape(bs.phi_face{1, 4} * U, [4 * bs.nfp, msh.nLElems]);

pl = computePressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = computePressure(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);
pd = computePressure(ud(1 : bs.nfp, :), ud(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ud(2 * bs.nfp + 1 : 3 * bs.nfp, :), ud(3 * bs.nfp + 1 : end, :), tc);
pu = computePressure(uu(1 : bs.nfp, :), uu(bs.nfp + 1 : 2 * bs.nfp, :), ...
    uu(2 * bs.nfp + 1 : 3 * bs.nfp, :), uu(3 * bs.nfp + 1 : end, :), tc);

p_min = min([p_cell; pl; pr; pd; pu], [], 1);
Cshift = max(0, (epsp - p_min) ./ (tc.gamma - 1));

U(1, 4 : 4 : end) = U(1, 4 : 4 : end) + Cshift;

if any(Cshift > 0)
    p_check = computeTestingPressure(U, tc, bs, msh);
    if any(~isfinite(p_check), 'all') || any(p_check < 0, 'all')
        error('Reviewer-recommended-method: pressure is negative or nonfinite after energy shift.');
    end
end

end

function p_all = computeTestingPressure(U,tc,bs,msh)
u_cell = bs.phi{1} * U;
p_cell = computePressure(u_cell(:, 1 : 4 : end), u_cell(:, 2 : 4 : end), ...
    u_cell(:, 3 : 4 : end), u_cell(:, 4 : 4 : end), tc);

ul = reshape(bs.phi_face{1, 1} * U, [4 * bs.nfp, msh.nLElems]);
ur = reshape(bs.phi_face{1, 2} * U, [4 * bs.nfp, msh.nLElems]);
ud = reshape(bs.phi_face{1, 3} * U, [4 * bs.nfp, msh.nLElems]);
uu = reshape(bs.phi_face{1, 4} * U, [4 * bs.nfp, msh.nLElems]);

pl = computePressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = computePressure(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);
pd = computePressure(ud(1 : bs.nfp, :), ud(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ud(2 * bs.nfp + 1 : 3 * bs.nfp, :), ud(3 * bs.nfp + 1 : end, :), tc);
pu = computePressure(uu(1 : bs.nfp, :), uu(bs.nfp + 1 : 2 * bs.nfp, :), ...
    uu(2 * bs.nfp + 1 : 3 * bs.nfp, :), uu(3 * bs.nfp + 1 : end, :), tc);

p_all = [p_cell; pl; pr; pd; pu];
end

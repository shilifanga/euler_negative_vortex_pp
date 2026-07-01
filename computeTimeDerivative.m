% Compute time derivative for the single-state P^k Euler formulation.
function Ut = computeTimeDerivative(msh, md, U, t, flux, tc, quad1, bs, IME)

% Evaluate physical quantities and flux functions on element faces.
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

Fl = computeF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl);
Fr = computeF(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pr);
Gd = computeG(ud(1 : bs.nfp, :), ud(bs.nfp + 1 : 2 * bs.nfp, :), ...
    ud(2 * bs.nfp + 1 : 3 * bs.nfp, :), ud(3 * bs.nfp + 1 : end, :), pd);
Gu = computeG(uu(1 : bs.nfp, :), uu(bs.nfp + 1 : 2 * bs.nfp, :), ...
    uu(2 * bs.nfp + 1 : 3 * bs.nfp, :), uu(3 * bs.nfp + 1 : end, :), pu);

% Element contributions.
if (bs.type == 200)
    Ut = zeros(1, 4 * msh.nLElems);
else
    u = bs.phi{1} * U;
    p = computePressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), ...
        u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);
    F = computeF(u(:, 1 : 4 : end), u(:, 2 : 4 : end), ...
        u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);
    G = computeG(u(:, 1 : 4 : end), u(:, 2 : 4 : end), ...
        u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);

    Ut = bs.phitw{2} * reshape(msh.elemJxix(:, msh.LElems) .* F, [bs.nep, 4 * msh.nLElems]) ...
        + bs.phitw{3} * reshape(msh.elemJetay(:, msh.LElems) .* G, [bs.nep, 4 * msh.nLElems]);
end

% Internal faces.
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

F_hat = computeInviscidFlux(ur(:, leLIDs), ul(:, reLIDs), pr(:, leLIDs), ...
    pl(:, reLIDs), Fr(:, leLIDs), Fl(:, reLIDs), 1, ones(1, nf), zeros(1, nf), flux, tc);
Ut(:, (-3 : 0)' + 4 * leLIDs) = Ut(:, (-3 : 0)' + 4 * leLIDs) ...
    - bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);
Ut(:, (-3 : 0)' + 4 * reLIDs) = Ut(:, (-3 : 0)' + 4 * reLIDs) ...
    + bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.intLFaces{2, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

F_hat = computeInviscidFlux(uu(:, leLIDs), ud(:, reLIDs), pu(:, leLIDs), ...
    pd(:, reLIDs), Gu(:, leLIDs), Gd(:, reLIDs), 2, zeros(1, nf), ones(1, nf), flux, tc);
Ut(:, (-3 : 0)' + 4 * leLIDs) = Ut(:, (-3 : 0)' + 4 * leLIDs) ...
    - bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);
Ut(:, (-3 : 0)' + 4 * reLIDs) = Ut(:, (-3 : 0)' + 4 * reLIDs) ...
    + bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

% Boundary faces.
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

yy = msh.elemCenter(2, leLIDs) + 0.5 * msh.elemLength(2, leLIDs) .* quad1.points;
u_ext = [tc.rho(msh.dm(1), yy, t); tc.m1(msh.dm(1), yy, t); ...
    tc.m2(msh.dm(1), yy, t); tc.E(msh.dm(1), yy, t)];
p_ext = computePressure(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), tc);
F_ext = computeF(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), p_ext);
F_hat = computeInviscidFlux(ul(:, leLIDs), u_ext, pl(:, leLIDs), p_ext, ...
    -Fl(:, leLIDs), -F_ext, 1, -ones(1, nf), zeros(1, nf), flux, tc);
Ut(:, (-3 : 0)' + 4 * leLIDs) = Ut(:, (-3 : 0)' + 4 * leLIDs) ...
    - bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{2, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

yy = msh.elemCenter(2, leLIDs) + 0.5 * msh.elemLength(2, leLIDs) .* quad1.points;
u_ext = [tc.rho(msh.dm(2), yy, t); tc.m1(msh.dm(2), yy, t); ...
    tc.m2(msh.dm(2), yy, t); tc.E(msh.dm(2), yy, t)];
p_ext = computePressure(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), tc);
F_ext = computeF(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), p_ext);
F_hat = computeInviscidFlux(ur(:, leLIDs), u_ext, pr(:, leLIDs), p_ext, ...
    Fr(:, leLIDs), F_ext, 1, ones(1, nf), zeros(1, nf), flux, tc);
Ut(:, (-3 : 0)' + 4 * leLIDs) = Ut(:, (-3 : 0)' + 4 * leLIDs) ...
    - bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{3, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

xx = msh.elemCenter(1, leLIDs) + 0.5 * msh.elemLength(1, leLIDs) .* quad1.points;
u_ext = [tc.rho(xx, tc.dm(3), t); tc.m1(xx, tc.dm(3), t); ...
    tc.m2(xx, tc.dm(3), t); tc.E(xx, tc.dm(3), t)];
p_ext = computePressure(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), tc);
F_ext = computeG(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), p_ext);
F_hat = computeInviscidFlux(ud(:, leLIDs), u_ext, pd(:, leLIDs), p_ext, ...
    -Gd(:, leLIDs), -F_ext, 2, zeros(1, nf), -ones(1, nf), flux, tc);
Ut(:, (-3 : 0)' + 4 * leLIDs) = Ut(:, (-3 : 0)' + 4 * leLIDs) ...
    - bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{4, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

xx = msh.elemCenter(1, leLIDs) + 0.5 * msh.elemLength(1, leLIDs) .* quad1.points;
u_ext = [tc.rho(xx, tc.dm(4), t); tc.m1(xx, tc.dm(4), t); ...
    tc.m2(xx, tc.dm(4), t); tc.E(xx, tc.dm(4), t)];
p_ext = computePressure(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), tc);
F_ext = computeG(u_ext(1 : bs.nfp, :), u_ext(bs.nfp + 1 : 2 * bs.nfp, :), ...
    u_ext(2 * bs.nfp + 1 : 3 * bs.nfp, :), u_ext(3 * bs.nfp + 1 : end, :), p_ext);
F_hat = computeInviscidFlux(uu(:, leLIDs), u_ext, pu(:, leLIDs), p_ext, ...
    Gu(:, leLIDs), F_ext, 2, zeros(1, nf), ones(1, nf), flux, tc);
Ut(:, (-3 : 0)' + 4 * leLIDs) = Ut(:, (-3 : 0)' + 4 * leLIDs) ...
    - bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

% Take care of the mass matrix.
Ut = (IME * Ut) ./ repelem(msh.elemJac(:, msh.LElems), 1, 4);

end

function F_hat = computeInviscidFlux(UL, UR, pl, pr, FL, FR, dir, nx, ny, flux, tc)
if flux ~= 4
    error('Only flux = 4 (HLLC) is supported')
end
F_hat = computeHLLCFlux(UL, UR, pl, pr, FL, FR, nx, ny, tc);
end

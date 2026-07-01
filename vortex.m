% Use EXRK-DG method to compute the two-dimensional Euler equations.
% This version uses one P^k polynomial state U only:
%   U(:,1:4:end) = rho, U(:,2:4:end) = m1,
%   U(:,3:4:end) = m2,  U(:,4:4:end) = E.

function L1err = vortex(N, tp, basisType, flux, pat, testcase, pp, postprocessPPType)

mtol = 1.0e-10;

if (nargin < 2)
    error('Not enough arguments')
end

if (nargin < 3) || isempty(basisType)
    basisType = 201;
end
if all(basisType ~= [201, 202])
    error('Only P1/P2 basis types 201 and 202 are supported')
end

if (nargin < 4) || isempty(flux)
    flux = 4;
end
if flux ~= 4
    error('Only flux = 4 (HLLC) is supported')
end

if (nargin < 5) || isempty(pat)
    pat = 0;
end
if (pat ~= 0 && pat ~= 1)
    error('Wrong pat index')
end

if (nargin < 8) || isempty(postprocessPPType)
    postprocessPPType = 'Standard-PP';
end

tc = createTestCase(1.4, testcase);

% P1/P2 use 3x3 or 4x4 Gauss--Lobatto element quadrature.
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 2, 102);
quad2 = GaussQuadratureRule_square([k + 2, k + 2], 202);

% basis function
bs = setBasisFunctionSet_square(quad1, quad2, basisType);

% Mesh.
msh = setRectMesh_rect(tc.dm, N, [2, 2, 2, 2], 201, 0);
md = computeMeshData_rect(msh);

% Inverse mass matrix on the reference square.
IME = inv(computeElementMatrix_refSquare([0, 0], [0, 0], bs, 1.0e-12, quad2));
IME = mychop(IME, mtol);

% Initial numerical solution.
U0 = computeInitialSolution_rect(msh, {tc.rho0, tc.m10, tc.m20, tc.E0}, bs, 1);
U0 = applypost2DEuler(U0, md, msh, bs, 0, tc, 0, testcase, pp, postprocessPPType);

if (pat == 1)
    plotDensity(msh, U0, quad2, bs);
end

cfl = setCFLNumber(k);

bt = getBTEXRK(k + 1);

t = 0;
while(t < tp - 1.0e-12)
    dt = setdt(msh, U0, t, tc, quad2, bs, cfl, tp);
    U0 = computeOneTimeStepEXRK(msh, md, U0, t, dt, flux, tc, ...
        quad1, bs, IME, bt, testcase, pp, postprocessPPType);
    t = t + dt;
 
end

L1err = computeError(msh, U0, t, tc, quad2, bs);

end

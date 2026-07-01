function U = computeOneTimeStepEXRK(msh, md, U0, t, dt, flux, tc, quad1, bs, IME, bt, testcase, pp, postprocessPPType)
%Time: SSPRK3
bs.deg = 2;
bt.c(2) = 1; bt.c(3) = 0.5;

if (nargin < 14) || isempty(postprocessPPType)
    postprocessPPType = 'Standard-PP';
end

switch bs.deg

    case 2
        LU = computeTimeDerivative(msh, md, U0, t, flux, tc, quad1, bs, IME);
        U = U0 + dt * LU;
        U = applypost2DEuler(U, md, msh, bs, dt, tc, t, testcase, pp, postprocessPPType);

        LU = computeTimeDerivative(msh, md, U, t + bt.c(2) * dt, flux, tc, quad1, bs, IME);
        U = 3 / 4 * U0 + 1 / 4 * (U + dt * LU);
        U = applypost2DEuler(U, md, msh, bs, dt, tc, t, testcase, pp, postprocessPPType);

        LU = computeTimeDerivative(msh, md, U, t + bt.c(3) * dt, flux, tc, quad1, bs, IME);
        U = 1 / 3 * U0 + 2 / 3 * (U + dt * LU);
        U = applypost2DEuler(U, md, msh, bs, dt, tc, t, testcase, pp, postprocessPPType);


    otherwise
        error('unsupported degree of polynomial')
end

end

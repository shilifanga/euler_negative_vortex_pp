function L1err = computeError(msh, U, t, tc, quad2, bs)
% Return only the L1 errors of density and pressure.

L1err = zeros(2, 1);

% Exact solution on the element quadrature points.
ct  = msh.elemCenter(:, msh.LElems);
h   = msh.elemLength(:, msh.LElems);
gpx = ct(1, :) + h(1, :) / 2 .* quad2.points(:, 1);
gpy = ct(2, :) + h(2, :) / 2 .* quad2.points(:, 2);
exact_rho = tc.rho(gpx, gpy, t);
exact_p   = tc.p(gpx, gpy, t);

% Numerical solution on the same quadrature points.
numerical_rho = bs.phi{1} * U(:, 1 : 4 : end);
numerical_m1  = bs.phi{1} * U(:, 2 : 4 : end);
numerical_m2  = bs.phi{1} * U(:, 3 : 4 : end);
numerical_E   = bs.phi{1} * U(:, 4 : 4 : end);
numerical_p   = computePressure(numerical_rho, numerical_m1, numerical_m2, numerical_E, tc);

J = msh.elemJac(:, msh.LElems);
L1err(1) = sum(J .* (quad2.weights' * abs(exact_rho - numerical_rho)));
L1err(2) = sum(J .* (quad2.weights' * abs(exact_p - numerical_p)));

end

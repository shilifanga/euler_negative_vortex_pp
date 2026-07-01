function tc = createTestCase(gamma,testcase)

switch testcase
    case 'normal-vortex'
        eps    = 5;
        % domain of computation
        tc.dm = [0, 10, 0, 10];

    case 'negative-vortex'
        eps    = 10.0828;%neg pressure
        % domain of computation
        tc.dm = [0, 10, 0, 10];

end
% ratio of specific heats, with gamma = 1.4 for air
tc.gamma = gamma;


% initial solution
rho_fs = 1;
u_fs   = 1;
v_fs   = 1;
p_fs   = 1;
T_fs   = p_fs / rho_fs;
x0     = 5;
y0     = 5;

du = @(x, y)eps / (2 * pi) * exp(0.5 * (1 - (x - x0).^2 - (y - y0).^2)) .* (y0 - y);
dv = @(x, y)eps / (2 * pi) * exp(0.5 * (1 - (x - x0).^2 - (y - y0).^2)) .* (x - x0);
dT = @(x, y)(1 - tc.gamma) * eps^2 / (8 * tc.gamma * pi^2) * exp(1 - (x - x0).^2 - (y - y0).^2);

tc.T0   = @(x, y)T_fs + dT(x, y);
tc.u0   = @(x, y)u_fs + du(x, y);
tc.v0   = @(x, y)v_fs + dv(x, y);
tc.rho0 = @(x, y)tc.T0(x, y).^(1 / (tc.gamma - 1));
tc.m10  = @(x, y)tc.rho0(x, y) .* tc.u0(x, y);
tc.m20  = @(x, y)tc.rho0(x, y) .* tc.v0(x, y);
tc.p0   = @(x, y)tc.T0(x, y).^(tc.gamma / (tc.gamma - 1));
tc.E0   = @(x, y)tc.p0(x, y) ./ (tc.gamma - 1) + 0.5 * tc.rho0(x, y) .* (tc.u0(x, y).^2 + tc.v0(x, y).^2);

% exact solution
tc.rho = @(x, y, t) tc.rho0(x - t, y - t);
tc.u   = @(x, y, t) tc.u0(x - t, y - t);
tc.v   = @(x, y, t) tc.v0(x - t, y - t);
tc.p   = @(x, y, t) tc.p0(x - t, y - t);
tc.m1  = @(x, y, t) tc.m10(x - t, y - t);
tc.m2  = @(x, y, t) tc.m20(x - t, y - t);
tc.E   = @(x, y, t) tc.E0(x - t, y - t);

 
end


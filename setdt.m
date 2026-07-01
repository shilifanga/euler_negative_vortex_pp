function dt = setdt(msh, U, t, tc, quad2, bs, cfl, tp)

% Evaluate means of density, momentum, energy, pressure and sound speed.
um = 0.25 * reshape(quad2.weights' * (bs.phi{1} * U), [4, msh.nLElems]);
pm = computePressure(um(1, :), um(2, :), um(3, :), um(4, :), tc);
cm = computeSpeedOfSound(um(1, :), pm, tc);

% Dissipation coefficient.
alpha = max(abs(um(2, :) ./ um(1, :)) + cm);
beta  = max(abs(um(3, :) ./ um(1, :)) + cm);

% Set the time step.
h  = msh.elemLength(:, msh.LElems);
dt = cfl / max(((alpha - 1) * (alpha > 1.0e-9) + 1) ./ h(1, :) + ((beta - 1) * (beta > 1.0e-9) + 1) ./ h(2, :));
if (t + dt) > tp
    dt = tp - t;
end

end

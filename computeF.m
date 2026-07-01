% the flux function in x direction 
function F = computeF(rho, m1, m2, E, p)

F = [m1; m1.^2 ./ rho + p; m1 .* m2 ./ rho; m1 ./ rho .* (E + p)];

end


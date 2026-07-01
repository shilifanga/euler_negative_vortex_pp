% equation of state
function p = computePressure(rho, m1, m2, E, tc)

p = (tc.gamma - 1) * (E - 0.5 * (m1.^2 + m2.^2) ./ rho);

 
end


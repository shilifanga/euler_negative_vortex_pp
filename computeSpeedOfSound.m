% speed of sound
function c = computeSpeedOfSound(rho, p, tc)

c = sqrt(tc.gamma * p ./ rho);

end

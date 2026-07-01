function H = computeHLLCFlux(UL, UR, pl, pr, FL, FR, nx, ny, tc)

nfp = size(UL, 1) / 4; % number of face points on each face
nf = size(UL, 2); % number of faces

% normal velocity
ul = nx .* (UL(nfp + 1 : 2 * nfp, :) ./ UL(1 : nfp, :)) + ny .* (UL(2 * nfp + 1 : 3 * nfp, :) ./ UL(1 : nfp, :));
ur = nx .* (UR(nfp + 1 : 2 * nfp, :) ./ UR(1 : nfp, :)) + ny .* (UR(2 * nfp + 1 : 3 * nfp, :) ./ UR(1 : nfp, :));

% speed of sound
al = computeSpeedOfSound(UL(1 : nfp, :), pl, tc);
ar = computeSpeedOfSound(UR(1 : nfp, :), pr, tc);

% wave speed
Sl = min(ul - al, ur - ar);
Sr = max(ul + al, ur + ar);
Sm = (UR(1 : nfp, :) .* ur .* (Sr - ur) - UL(1 : nfp, :) .* ul .* (Sl - ul) + pl - pr) ./ (UR(1 : nfp, :) .* (Sr - ur) - UL(1 : nfp, :) .* (Sl - ul));

% intermediate pressure
ps = UL(1 : nfp, :) .* (Sl - ul) .* (Sm - ul) + pl;

% other quantities
cl = 0.5 * (1 + (abs(Sl) - abs(Sm)) ./ (Sl - Sm)); 
cr = 0.5 * (1 - (abs(Sr) - abs(Sm)) ./ (Sr - Sm));
vp = [zeros(nfp, nf); repelem(nx, nfp, 1); repelem(ny, nfp, 1); Sm];

% HLLC flux
H = repmat(cl, [4, 1]) .* FL + repmat(cr, [4, 1]) .* FR - repmat((cl + cr - 1) .* ps, [4, 1]) .* vp + repmat(0.5 * (abs(Sl) - Sl .* (abs(Sl) - abs(Sm)) ./ (Sl - Sm)), [4, 1]) .* UL - repmat(0.5 * (abs(Sr) - Sr .* (abs(Sr) - abs(Sm)) ./ (Sr - Sm)), [4, 1]) .* UR;

end
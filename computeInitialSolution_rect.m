% msh    : mesh of rectangular element in 2D
% u0     : cell array of exact initial solutions
% bs     : basis function set data or degree 1/2
% layout : 0 or 1, layout of U0

function U0 = computeInitialSolution_rect(msh, u0, bs, layout)

if (nargin < 2)
    error('Not enough arguments')
end
if (msh.type ~= 201) && (msh.type ~= 202)
    error('Wrong mesh type')
end
if (nargin < 3) || isempty(bs)
    bs = 1;
end

if isfloat(bs) && isscalar(bs)
    k = floor(bs);
    if all(k ~= [1, 2])
        error('Only P1/P2 basis degrees 1 and 2 are supported')
    end
    quad1 = GaussQuadratureRule_line(k + 2, 102);
    quad2 = GaussQuadratureRule_square([k + 2, k + 2], 202);
    bs = setBasisFunctionSet_square(quad1, quad2, 200 + k);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'square')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    quad2 = GaussQuadratureRule_square(bs.neps, bs.elemPointsType);
else
    error('Wrong argument bs')
end

if (nargin < 4) || isempty(layout)
    layout = 0;
end
if (layout ~= 0) && (layout ~= 1)
    error('Wrong argument layout')
end

ME = computeElementMatrix_refSquare([0, 0], [0, 0], bs, 1.0e-12, quad2);

ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);

nv = length(u0);
U0 = zeros(bs.nb, nv * msh.nLElems);
switch layout
    case 0
        for i = 1 : nv
            cols = (i - 1) * msh.nLElems + 1 : i * msh.nLElems;
            U0(:, cols) = ME \ (bs.phitw{1} * u0{i}( ...
                ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1), ...
                ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2)));
        end
    case 1
        for i = 1 : nv
            U0(:, i : nv : end) = ME \ (bs.phitw{1} * u0{i}( ...
                ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1), ...
                ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2)));
        end
end

end
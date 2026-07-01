% faceOrderType = 1 : the four faces are ordered by left, right, bottom and
% top, and the points on the four face are placed in an increasing order
% w.r.t. x or y.
% faceOrderType = 2 : the four faces are ordered by left, right, bottom and
% top, and the points on the four face are placed in a counterclockwise
% order.

function bs = setBasisFunctionSet_square(quad1, quad2, basisType, faceOrderType)

if (nargin < 3)
    error('Not enough arguments')
end
if all(basisType ~= [201, 202])
    error('Only P1/P2 basis types 201 and 202 are supported')
end

if (nargin < 4) || isempty(faceOrderType)
    faceOrderType = 1;
end
if (faceOrderType ~= 1) && (faceOrderType ~= 2)
    error('Wrong face points order')
end

bs.refGeom = 'square';
bs.type = basisType;
bs.deg = mod(basisType, 10);
bs.elemPointsType = quad2.type;
bs.facePointsType = quad1.type;
bs.faceOrderType = faceOrderType;

% Basis functions and first derivatives at element quadrature nodes.
bs.phi = {basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), basisType, [0, 0]), ...
          basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), basisType, [1, 0]), ...
          basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), basisType, [0, 1])};

% Basis functions at face quadrature nodes.
bs.phi_face = cell(2, 4);
if (faceOrderType == 1)
    bs.phi_face(1, :) = {basisFunctionSet_square(-1, quad1.points, basisType, [0, 0]), ...
                         basisFunctionSet_square( 1, quad1.points, basisType, [0, 0]), ...
                         basisFunctionSet_square(quad1.points, -1, basisType, [0, 0]), ...
                         basisFunctionSet_square(quad1.points,  1, basisType, [0, 0])};
else
    bs.phi_face(1, :) = {basisFunctionSet_square(-1, flip(quad1.points), basisType, [0, 0]), ...
                         basisFunctionSet_square( 1, quad1.points,       basisType, [0, 0]), ...
                         basisFunctionSet_square(quad1.points,      -1,  basisType, [0, 0]), ...
                         basisFunctionSet_square(flip(quad1.points), 1,  basisType, [0, 0])};
end
bs.phi_face(2, :) = {flip(bs.phi_face{1, 1}), flip(bs.phi_face{1, 2}), ...
                     flip(bs.phi_face{1, 3}), flip(bs.phi_face{1, 4})};

% Weighted basis data for element and face integration.
bs.phitw = {(quad2.weights .* bs.phi{1})', ...
            (quad2.weights .* bs.phi{2})', ...
            (quad2.weights .* bs.phi{3})'};

bs.phitw_face = cell(2, 4);
for i = 1 : 2
    for j = 1 : 4
       bs.phitw_face{i, j} = (quad1.weights .* bs.phi_face{i, j})';
    end
end

bs.nb = size(bs.phi{1}, 2);
bs.nb2 = bs.nb * bs.nb;
bs.nep = quad2.np;
bs.nfp = quad1.np;
bs.neps = quad2.nps;

end
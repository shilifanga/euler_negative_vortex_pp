function quad = GaussQuadratureRule_square(np, type)
% Tensor-product Gauss--Lobatto quadrature with at most 4 x 4 points.

if (numel(np) ~= 2)
    error('The first argument should be a 2D vector')
end
if type ~= 202
    error('Only tensor-product Gauss--Lobatto quadrature type 202 is supported')
end
if any(np < 2) || any(np > 4)
    error('Only 2x2, 3x3, or 4x4 Gauss--Lobatto quadrature is supported')
end

quad.type = type;
quad.np = prod(np);
quad.nps = np;

quad1 = GaussQuadratureRule_line(np(1), 102);
quad2 = GaussQuadratureRule_line(np(2), 102);
quad.points = [repmat(quad1.points, [np(2), 1]), repelem(quad2.points, np(1))];
quad.weights = reshape(quad1.weights * quad2.weights', [quad.np, 1]);

end
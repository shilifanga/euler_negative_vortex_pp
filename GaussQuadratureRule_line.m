function quad = GaussQuadratureRule_line(np, type)
% Gauss--Lobatto quadrature with at most four points.

if type ~= 102
    error('Only Gauss--Lobatto quadrature type 102 is supported')
end

quad.type = type;
quad.np = np;

switch np
    case 2
        quad.points = [-1; 1];
        quad.weights = [1; 1];
    case 3
        quad.points = [-1; 0; 1];
        quad.weights = [1 / 3; 4 / 3; 1 / 3];
    case 4
        quad.points = [-1; -sqrt(1 / 5); sqrt(1 / 5); 1];
        quad.weights = [1 / 6; 5 / 6; 5 / 6; 1 / 6];
    otherwise
        error('Only 2, 3, or 4 Gauss--Lobatto points are supported')
end

end
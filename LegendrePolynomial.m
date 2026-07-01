function r = LegendrePolynomial(x, poly_deg, der_order)
% Orthogonal modal polynomials used by the P1/P2 basis.

switch der_order
    case 0
        switch poly_deg
            case 0
                r = ones(size(x));
            case 1
                r = x;
            case 2
                r = x.^2 - 1 / 3;
            otherwise
                error('Only polynomial degrees 0, 1, and 2 are supported')
        end
    case 1
        switch poly_deg
            case 0
                r = zeros(size(x));
            case 1
                r = ones(size(x));
            case 2
                r = 2 * x;
            otherwise
                error('Only polynomial degrees 0, 1, and 2 are supported')
        end
    otherwise
        error('Only derivative orders 0 and 1 are supported')
end

end
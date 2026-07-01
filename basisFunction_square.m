function r = basisFunction_square(x, y, type, index, der_order)
% Orthogonal P1/P2 basis functions on the reference square.

if numel(der_order) ~= 2
    error('The argument der_order should be a 2D vector')
end

switch type
    case {201, 202}
        [i, j] = getSplitIndex(type, index);
        r = LegendrePolynomial(x, i, der_order(1)) .* LegendrePolynomial(y, j, der_order(2));
    otherwise
        error('Only P1/P2 basis types 201 and 202 are supported')
end

    function [i, j] = getSplitIndex(localType, localIndex)
        switch localType
            case 201
                switch localIndex
                    case 1
                        i = 0; j = 0;
                    case 2
                        i = 1; j = 0;
                    case 3
                        i = 0; j = 1;
                    otherwise
                        error('P1 has three basis functions')
                end
            case 202
                switch localIndex
                    case 1
                        i = 0; j = 0;
                    case 2
                        i = 1; j = 0;
                    case 3
                        i = 0; j = 1;
                    case 4
                        i = 2; j = 0;
                    case 5
                        i = 1; j = 1;
                    case 6
                        i = 0; j = 2;
                    otherwise
                        error('P2 has six basis functions')
                end
        end
    end

end
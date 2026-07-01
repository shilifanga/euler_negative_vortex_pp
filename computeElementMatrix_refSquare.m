% trial_der : derivative order of trial basis functions
% test_der  : derivative order of test basis functions
% bs        : basis function set data or degree 1/2
% mtol      : tolerance for the matrix entries

function EM = computeElementMatrix_refSquare(trial_der, test_der, bs, mtol, quad2)

if (nargin < 1) || isempty(trial_der)
    trial_der = [0, 0];
end
if (nargin < 2) || isempty(test_der)
    test_der = [0, 0];
end
if (length(trial_der) ~= 2) || (length(test_der) ~= 2)
    error('Wrong size of derivative order')
end
if (nargin < 3) || isempty(bs)
    bs = 1;
end
if (nargin < 4) || isempty(mtol)
    mtol = 1.0e-12;
end

check_derivative(trial_der);
check_derivative(test_der);

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
    if (nargin < 5) || isempty(quad2)
        quad2 = GaussQuadratureRule_square(bs.neps, bs.elemPointsType);
    end
else
    error('Wrong argument bs')
end

phi_trial = basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), bs.type, trial_der);
phi_test = basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), bs.type, test_der);
phitw_test = (quad2.weights .* phi_test)';
EM = phitw_test * phi_trial;
EM = mychop(EM, mtol);

end

function check_derivative(der)
if all(der(:) == [0; 0]) || all(der(:) == [1; 0]) || all(der(:) == [0; 1])
    return
end
error('Only derivative orders [0,0], [1,0], and [0,1] are supported')
end
% pp = 0 : no limiter
% pp = 1 : pp limiter
  
function U = applyPPLimiter3(U, pp, tc, bs1, limiterMode)

if (nargin < 2) || isempty(pp)
    pp = 1;
end
if (nargin < 5) || isempty(limiterMode)
    limiterMode = 'Standard-PP';
end

if (pp == 0) || (bs1.type == 200)
    return
end

epsrho = 1.0e-12;
epsp = 1.0e-12;

% Step 1. Density positivity. Only density high-order moments are changed.
rho_bar = U(1, 1 : 4 : end);
rho_pts = bs1.phi{1} * U(:, 1 : 4 : end);
rho_min = min(rho_pts, [], 1);

theta_rho = ones(1, size(rho_pts, 2));
viol = rho_min < epsrho;
den = rho_bar(viol) - rho_min(viol);
den = max(den, 1.0e-300);
theta_rho(viol) = (rho_bar(viol) - epsrho) ./ den;
theta_rho = min(1, max(0, theta_rho));

U(2 : end, 1 : 4 : end) = theta_rho .* U(2 : end, 1 : 4 : end);

rhoPP = bs1.phi{1} * U(:, 1 : 4 : end);
rhoPP_min = min(rhoPP, [], 1);
bad_rho = rhoPP_min < epsrho;
if any(bad_rho)
    den = rho_bar(bad_rho) - rhoPP_min(bad_rho);
    den = max(den, 1.0e-300);
    theta_fix = (rho_bar(bad_rho) - epsrho) ./ den;
    theta_fix = min(1, max(0, theta_fix));
    bad_cols = find(bad_rho);
    U(2 : end, 4 * bad_cols - 3) = theta_fix .* U(2 : end, 4 * bad_cols - 3);
    rhoPP = bs1.phi{1} * U(:, 1 : 4 : end);
end
if any(rhoPP < 0, 'all')
    error('applyPPLimiter3: density is still negative after limiting.');
end

limiterPP = strcmp(limiterMode, 'Standard-PP');
if limiterPP
    % Step 2. Pressure positivity. By default use the standard Zhang-Shu
    % quadratic root. The older linear scaling remains as a switch option.
    pressureLimiterType = getPressureLimiterType();

    u = bs1.phi{1} * U;
    p = computePressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), ...
        u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);

    pbar = computePressure(U(1, 1 : 4 : end), U(1, 2 : 4 : end), ...
        U(1, 3 : 4 : end), U(1, 4 : 4 : end), tc);

    if any(pbar <= epsp)
        error('applyPPLimiter3 needs positive cell-average pressure.');
    end

    switch pressureLimiterType
        case 'Zhang-ShuQuadratic'
            t_alpha = ones(size(p));
            pos = find(p <= epsp);
            if ~isempty(pos)
                col = ceil(pos / size(p, 1));
                rhom_pos = U(1, 4 * col - 3)';
                m1m_pos  = U(1, 4 * col - 2)';
                m2m_pos  = U(1, 4 * col - 1)';
                Em_pos   = U(1, 4 * col)';
                pm_pos   = pbar(col)';
                rho_pos  = getMatEntries(u(:, 1 : 4 : end), pos);
                m1_pos   = getMatEntries(u(:, 2 : 4 : end), pos);
                m2_pos   = getMatEntries(u(:, 3 : 4 : end), pos);
                E_pos    = getMatEntries(u(:, 4 : 4 : end), pos);
                p_pos    = p(pos);

                AA = (tc.gamma - 1) * (m1m_pos .* m1_pos + m2m_pos .* m2_pos ...
                    - rhom_pos .* E_pos - rho_pos .* Em_pos);
                aa = rhom_pos .* pm_pos + rho_pos .* p_pos + AA;
                bb = -AA - 2 * rhom_pos .* pm_pos - epsp * (rho_pos - rhom_pos);
                cc = rhom_pos .* (pm_pos - epsp);

                disc = max(bb .^ 2 - 4 * aa .* cc, 0);
                sqrt_disc = sqrt(disc);
                root_direct = (-bb - sqrt_disc) ./ (2 * aa);
                root_conj = 2 * cc ./ (-bb + sqrt_disc);
                root_linear = -cc ./ bb;

                valid_direct = isfinite(root_direct) & (root_direct >= 0) & (root_direct <= 1);
                valid_conj = isfinite(root_conj) & (root_conj >= 0) & (root_conj <= 1);
                valid_linear = isfinite(root_linear) & (root_linear >= 0) & (root_linear <= 1);

                t_root = root_direct;
                both_valid = valid_direct & valid_conj;
                t_root(both_valid) = min(root_direct(both_valid), root_conj(both_valid));
                t_root(~valid_direct & valid_conj) = root_conj(~valid_direct & valid_conj);
                t_root(~valid_direct & ~valid_conj & valid_linear) = ...
                    root_linear(~valid_direct & ~valid_conj & valid_linear);
                t_root(~isfinite(t_root)) = 0;
                t_alpha(pos) = min(1, max(0, t_root));
            end
            theta_p = min(t_alpha, [], 1);

        case 'ConcavityLinear'
            theta_p_pts = ones(size(p));
            bad = p <= epsp;
            pbar_pts = repmat(pbar, size(p, 1), 1);
            theta_p_pts(bad) = (pbar_pts(bad) - epsp) ./ ...
                (pbar_pts(bad) - p(bad) + epsp);
            theta_p = min(theta_p_pts, [], 1);

        otherwise
            error('Unknown pressureLimiterType')
    end
    theta_p = min(1, max(0, theta_p));

    % Shrink all high-order conservative moments from the cell mean.
    U(2 : end, :) = repelem(theta_p, 1, 4) .* U(2 : end, :);

    u = bs1.phi{1} * U;
    p_after = computePressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), ...
        u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);

    if any(p_after < 0, 'all')
        error('applyPPLimiter3: pressure is still negative after limiting.');
    end
end

end

function entries = getMatEntries(mat, pos)
entries = mat(pos);
end

function pressureLimiterType = getPressureLimiterType()
pressureLimiterType = getenv('EULER_PP_PRESSURE_LIMITER');
if isempty(pressureLimiterType)
    pressureLimiterType = 'Zhang-ShuQuadratic';
end
end

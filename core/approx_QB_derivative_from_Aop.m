function [dQB, QB_left, QB_right, info] = approx_QB_derivative_from_Aop(A_op, mB, lam0, h, opts)
% Approximate d/dlambda QB(lambda) using a shared selected column set.
%
% Workflow:
% 1. build A(lambda-h), A(lambda+h)
% 2. do pivoted QR at both points and truncate
% 3. intersect the selected original columns
% 4. rebuild QB on that shared column set
% 5. return centered finite-difference derivative

    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    lam_left = lam0 - h;
    lam_right = lam0 + h;

    A_left = local_maybe_normalize(A_op(lam_left), opts);
    A_right = local_maybe_normalize(A_op(lam_right), opts);

    [cols_left, rank_left, gap_left] = local_select_columns(A_left, opts);
    [cols_right, rank_right, gap_right] = local_select_columns(A_right, opts);

    common_cols = intersect(cols_left, cols_right, 'stable');

    if isempty(common_cols)
        QB_left = zeros(mB, 0, class(A_left));
        QB_right = zeros(mB, 0, class(A_right));
        dQB = zeros(mB, 0, class(A_left));
    else
        QB_left = local_build_QB(A_left(:, common_cols), mB, opts);
        QB_right = local_build_QB(A_right(:, common_cols), mB, opts);
        dQB = (QB_right - QB_left) / (2*h);
    end

    if nargout >= 4
        info = struct();
        info.rank_left = rank_left;
        info.rank_right = rank_right;
        info.n_common = numel(common_cols);
        info.gap_left = gap_left;
        info.gap_right = gap_right;
        info.common_cols = common_cols;
    end
end


function [selected_cols, rank_qr, gap] = local_select_columns(A, opts)
    use_pivot = local_get_opt(opts, 'pivot', true);
    qr_tau = local_get_opt(opts, 'qr_tau', []);

    if use_pivot
        [~, R, perm] = local_qr_pivot(A);
    else
        [~, R] = qr(A, 0);
        perm = 1:size(A, 2);
    end

    diagR = abs(diag(R));
    rank_qr = local_rank_from_diag(diagR, qr_tau);
    selected_cols = perm(1:rank_qr);
    gap = local_gap(diagR, rank_qr);
end


function QB = local_build_QB(A, mB, opts)
    [Q, R] = qr(A, 0);

    if local_get_opt(opts, 'sign_fix', true)
        d = sign(diag(R));
        d(d == 0) = 1;
        Q = Q * diag(d);
    end

    QB = Q(1:mB, :);
end


function A = local_maybe_normalize(A, opts)
    if local_get_opt(opts, 'normalize_columns', false)
        A = A ./ max(vecnorm(A), 1e-300);
    end
end


function [Q, R, perm] = local_qr_pivot(A)
    try
        [Q, R, perm] = qr(A, 0, 'vector');
    catch
        [Q, R, E] = qr(A, 0);
        perm = local_permvec(E);
    end
end


function perm = local_permvec(E)
    if isvector(E)
        perm = E(:).';
        return;
    end

    [~, perm] = max(abs(E), [], 1);
    perm = perm(:).';
end


function rank_qr = local_rank_from_diag(diagR, qr_tau)
    if isempty(diagR)
        rank_qr = 0;
        return;
    end

    if isempty(qr_tau)
        rank_qr = numel(diagR);
        return;
    end

    rank_qr = find(diagR >= qr_tau * diagR(1), 1, 'last');
    if isempty(rank_qr)
        rank_qr = 0;
    end
end


function gap = local_gap(diagR, rank_qr)
    if isempty(diagR) || rank_qr < 1 || rank_qr >= numel(diagR)
        gap = NaN;
        return;
    end

    gap = diagR(rank_qr) / max(diagR(rank_qr + 1), eps(class(diagR)));
end


function value = local_get_opt(opts, name, default_value)
    if isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    else
        value = default_value;
    end
end

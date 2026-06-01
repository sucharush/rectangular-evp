function [QB, info, R] = build_QB_from_A(A, mB, opts)
% Build the boundary block QB from a full collocation matrix A.
%
% opts.normalize_columns : true/false
% opts.pivot             : true/false
% opts.sign_fix          : true/false
% opts.qr_tau            : truncation threshold on abs(diag(R))

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    normalize_columns = get_logical_opt(opts, 'normalize_columns', false);
    pivot = get_logical_opt(opts, 'pivot', false);
    sign_fix = get_logical_opt(opts, 'sign_fix', false);
    qr_tau = get_opt(opts, 'qr_tau', []);

    if normalize_columns
        column_norms = max(vecnorm(A), 1e-300);
        A = A ./ column_norms;
    else
        column_norms = [];
    end

    if pivot
        [Q, R, perm] = local_qr_pivot(A);
    else
        ncol = size(A, 2);
        [Q, R] = qr(A, 0);
        perm = 1:ncol;
    end

    if sign_fix
        d = sign(diag(R));
        d(d == 0) = 1;
        D = diag(d);
        Q = Q * D;
        R = D * R;
    end
    % condR = cond(R);

    diagR = abs(diag(R));
    rank_qr = local_truncation_rank(diagR, qr_tau);

    Q = Q(:, 1:rank_qr);
    QB = Q(1:mB, :);
    Rp = R(1:rank_qr, 1:rank_qr);
    condR = cond(Rp);


    info = struct();
    info.perm = perm;
    info.diagR = diagR;
    info.rank = rank_qr;
    info.selected_cols = perm(1:rank_qr);
    info.column_norms = column_norms;
end


function perm = local_permvec(E)
    if isvector(E)
        perm = E(:).';
        return;
    end

    [~, perm] = max(abs(E), [], 1);
    perm = perm(:).';
end


function [Q, R, perm] = local_qr_pivot(A)
    try
        [Q, R, perm] = qr(A, 0, 'vector');
    catch
        [Q, R, E] = qr(A, 0);
        perm = local_permvec(E);
    end
end


function rank_qr = local_truncation_rank(diagR, qr_tau)
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


function value = get_opt(opts, name, default_value)
    if isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    else
        value = default_value;
    end
end


function value = get_logical_opt(opts, name, default_value)
    value = get_opt(opts, name, default_value);
    value = logical(value);
end

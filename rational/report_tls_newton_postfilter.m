function report_tls_newton_postfilter(out, pencil_data, lam_true, top_k, label)
% Block-based reduced-subspace post-filter for the Newton-form linearized pencil.
% Reconstructs each TLS eigenvector's block subspace via the Newton basis prefix,
% then scores candidates by the smallest singular value of R(lambda) restricted to
% that subspace.

    tls_data = pencil_data.tls;
    lam_tls_all = tls_data.lam_all;
    vec_tls_all = tls_data.vec_all;
    idx_keep = tls_data.idx_keep;
    m = out.m;
    q = out.q;

    if size(vec_tls_all, 1) ~= m * q
        error('report_tls_newton_postfilter: TLS eigenvector length does not match m*q.');
    end

    nkeep = numel(idx_keep);
    reduced_score = nan(nkeep, 1);
    first_block_score = nan(nkeep, 1);
    matched_true = nan(nkeep, 1);
    err_to_true = nan(nkeep, 1);
    subspace_dim = nan(nkeep, 1);
    n_valid_blocks = nan(nkeep, 1);
    rel_pencil_res = tls_data.rel_res_all(idx_keep);
    lam_keep = lam_tls_all(idx_keep);

    for t = 1:nkeep
        j = idx_keep(t);
        lam = lam_keep(t);
        x = vec_tls_all(:, j);

        bj = local_newton_basis_prefix(lam, out.sigma, out.beta, out.h, out.k);
        [U, n_valid] = local_build_tls_block_subspace(x, bj, q);
        if isempty(U)
            continue;
        end

        Rb = aaa_eval_matrix_barycentric(lam, out.zj, out.wj, out.D);
        RU = Rb * U;
        svals_reduced = svd(RU, 'econ');
        if isempty(svals_reduced)
            continue;
        end

        reduced_score(t) = svals_reduced(end);
        subspace_dim(t) = size(U, 2);
        n_valid_blocks(t) = n_valid;

        y_first = local_extract_first_block_y(x, q);
        if ~isempty(y_first)
            y_first_norm = norm(y_first);
            first_block_score(t) = norm(Rb * y_first) / max(y_first_norm, 1e-14);
        end

        [err_to_true(t), idx_true] = min(abs(lam - lam_true(:)));
        matched_true(t) = lam_true(idx_true);
    end

    valid_score_idx = find(isfinite(reduced_score));
    [~, order] = sort(reduced_score(valid_score_idx), 'ascend');
    idx_ranked = valid_score_idx(order);
    idx_top = idx_ranked(1:min(top_k, numel(idx_ranked)));

    fprintf('\nTLS reduced-subspace post-filter: %s\n', label);
    fprintf('  prefilter: real(lambda) in [%.6f, %.6f], |imag(lambda)| <= %.3g\n', ...
        tls_data.interval(1), tls_data.interval(2), tls_data.imag_tol);
    fprintf('  prefiltered %d of %d TLS pairs; reporting top %d by reduced score.\n', ...
        nkeep, numel(lam_tls_all), numel(idx_top));

    if isempty(idx_top)
        fprintf('  no prefiltered TLS pairs produced a valid reduced subspace score.\n');
    else
        fprintf('%4s  %16s  %11s  %12s  %12s  %12s  %6s  %6s  %16s  %12s\n', ...
            'rank', 'real(lambda)', 'imag(lambda)', 'red_score', 'first_score', ...
            'pencil_res', 'dimU', 'nblk', 'matched_true', 'err_true');

        for t = 1:numel(idx_top)
            j = idx_top(t);
            fprintf('%4d  %16.12f  %11.3e  %12.3e  %12.3e  %12.3e  %6d  %6d  %16.12f  %12.3e\n', ...
                t, real(lam_keep(j)), imag(lam_keep(j)), reduced_score(j), ...
                first_block_score(j), rel_pencil_res(j), round(subspace_dim(j)), ...
                round(n_valid_blocks(j)), matched_true(j), err_to_true(j));
        end
    end
end


function bj = local_newton_basis_prefix(z, sigma, beta, h, k)
    m = numel(beta);
    bj = ones(m, 1);

    for j = 1:(m-1)
        denom = beta(j) * (h(j) - k(j) * z);
        if ~isfinite(denom) || abs(denom) < 1e-14
            bj(j+1:end) = NaN;
            return;
        end
        bj(j+1) = ((z - sigma(j)) / denom) * bj(j);
    end
end


function [U, n_valid] = local_build_tls_block_subspace(x, bj, q)
    m = numel(bj);
    X = reshape(x, q, m);
    Y = zeros(q, m);
    keep = false(1, m);

    for j = 1:m
        if ~isfinite(bj(j)) || abs(bj(j)) < 1e-14
            continue;
        end

        yj = X(:, j) / bj(j);
        if any(~isfinite(yj))
            continue;
        end

        Y(:, j) = yj;
        keep(j) = true;
    end

    n_valid = nnz(keep);
    if n_valid == 0
        U = [];
        return;
    end

    Y = Y(:, keep);
    [Q, R] = qr(Y, 0);
    diagR = abs(diag(R));

    if isempty(diagR)
        U = [];
        n_valid = 0;
        return;
    end

    tol = max(size(R)) * eps(max(diagR));
    rankU = nnz(diagR > tol);
    U = Q(:, 1:rankU);

    if isempty(U)
        n_valid = 0;
    end
end


function y = local_extract_first_block_y(x, q)
    if numel(x) < q
        y = [];
        return;
    end

    y = x(1:q);
    y_norm = norm(y);
    if ~isfinite(y_norm) || y_norm == 0
        y = [];
    end
end

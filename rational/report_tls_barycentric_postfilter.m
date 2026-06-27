function report_tls_barycentric_postfilter(out, pencil_data, lam_true, top_k, label)
% Block-based reduced-subspace post-filter for the star-barycentric linearized
% pencil. Reconstructs each TLS eigenvector's block subspace from the star nodes,
% then scores candidates by the smallest singular value of R(lambda) restricted to
% that subspace.

    tls_data = pencil_data.tls;
    lam_tls_all = tls_data.lam_all;
    vec_tls_all = tls_data.vec_all;
    idx_keep = tls_data.idx_keep;

    nblocks = numel(out.zj);
    q = out.q;

    if size(vec_tls_all, 1) ~= (nblocks + 1) * q
        error('report_tls_barycentric_postfilter: TLS eigenvector length does not match (numel(zj)+1)*q.');
    end

    nkeep = numel(idx_keep);
    reduced_score = nan(nkeep, 1);
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

        [U, n_valid] = local_build_star_barycentric_tls_subspace(x, lam, out.zj, q);
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
        fprintf('%4s  %16s  %11s  %12s  %12s  %6s  %6s  %16s  %12s\n', ...
            'rank', 'real(lambda)', 'imag(lambda)', 'red_score', ...
            'pencil_res', 'dimU', 'nblk', 'matched_true', 'err_true');

        for t = 1:numel(idx_top)
            j = idx_top(t);
            fprintf('%4d  %16.12f  %11.3e  %12.3e  %12.3e  %6d  %6d  %16.12f  %12.3e\n', ...
                t, real(lam_keep(j)), imag(lam_keep(j)), reduced_score(j), ...
                rel_pencil_res(j), round(subspace_dim(j)), ...
                round(n_valid_blocks(j)), matched_true(j), err_to_true(j));
        end
    end
end


function [U, n_valid] = local_build_star_barycentric_tls_subspace(x, lam, z_nodes, q)
    nblocks = numel(z_nodes);
    X = reshape(x, q, nblocks + 1);
    Y = zeros(q, nblocks + 1);
    keep = false(1, nblocks + 1);

    % The first block is the global x block.
    x0 = X(:, 1);
    if all(isfinite(x0)) && norm(x0) > 0
        Y(:, 1) = x0;
        keep(1) = true;
    end

    for j = 1:nblocks
        scale = lam - z_nodes(j);
        if ~isfinite(scale) || abs(scale) < 1e-14
            continue;
        end

        xj = scale * X(:, j + 1);
        if any(~isfinite(xj)) || norm(xj) == 0
            continue;
        end

        Y(:, j + 1) = xj;
        keep(j + 1) = true;
    end

    n_valid = nnz(keep);
    if n_valid == 0
        U = [];
        return;
    end

    Ykeep = Y(:, keep);
    [Q, R] = qr(Ykeep, 0);
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

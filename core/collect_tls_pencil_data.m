function tls_data = collect_tls_pencil_data(A, B, interval, imag_tol)
% Collect TLS eigenpairs of the rectangular pencil A - lambda B and prefilter
% them to a real interval with a small imaginary tolerance. Returns residuals
% and the kept index set for downstream ranking/reporting.

    [vec_tls_all, lam_tls_all] = tls_pencil_eigs(full(A), full(B));

    lam_tls_all = lam_tls_all(:);
    nlam = numel(lam_tls_all);

    if size(vec_tls_all, 2) ~= nlam
        error('collect_tls_pencil_data: inconsistent TLS eigenvector/eigenvalue sizes.');
    end

    abs_res = nan(nlam, 1);
    rel_res = nan(nlam, 1);

    normA = norm(A, 'fro');
    normB = norm(B, 'fro');

    for j = 1:nlam
        lam = lam_tls_all(j);
        x = vec_tls_all(:, j);

        if ~all(isfinite([real(lam), imag(lam)]))
            continue;
        end

        xnorm = norm(x);
        if xnorm == 0 || ~isfinite(xnorm)
            continue;
        end

        r = (A - lam * B) * x;
        abs_res(j) = norm(r);
        rel_res(j) = abs_res(j) / max((normA + abs(lam) * normB) * xnorm, 1e-14);
    end

    keep_mask = isfinite(real(lam_tls_all)) ...
        & isfinite(imag(lam_tls_all)) ...
        & isfinite(rel_res) ...
        & real(lam_tls_all) >= interval(1) ...
        & real(lam_tls_all) <= interval(2) ...
        & abs(imag(lam_tls_all)) <= imag_tol;

    tls_data = struct();
    tls_data.A = A;
    tls_data.B = B;
    tls_data.lam_all = lam_tls_all;
    tls_data.vec_all = vec_tls_all;
    tls_data.abs_res_all = abs_res;
    tls_data.rel_res_all = rel_res;
    tls_data.keep_mask = keep_mask;
    tls_data.idx_keep = find(keep_mask);
    tls_data.interval = interval;
    tls_data.imag_tol = imag_tol;
end

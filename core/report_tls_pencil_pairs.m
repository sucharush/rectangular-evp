function report_tls_pencil_pairs(tls_data, lam_true, top_k, label)
% Print the prefiltered TLS pencil candidates ranked by relative residual,
% matched to the nearest reference eigenvalue in lam_true.

    idx_keep = tls_data.idx_keep;
    rel_res = tls_data.rel_res_all;
    lam_tls_all = tls_data.lam_all;
    abs_res = tls_data.abs_res_all;

    [~, order] = sort(rel_res(idx_keep), 'ascend');
    idx_ranked = idx_keep(order);
    idx_top = idx_ranked(1:min(top_k, numel(idx_ranked)));

    fprintf('\nTLS pencil candidates: %s\n', label);
    fprintf('  prefilter: real(lambda) in [%.6f, %.6f], |imag(lambda)| <= %.3g\n', ...
        tls_data.interval(1), tls_data.interval(2), tls_data.imag_tol);
    fprintf('  kept %d of %d TLS pairs; reporting top %d by relative residual.\n', ...
        numel(idx_keep), numel(lam_tls_all), numel(idx_top));

    if isempty(idx_top)
        fprintf('  no TLS pairs passed the prefilter.\n');
    else
        fprintf('%4s  %16s  %11s  %12s  %12s  %16s  %12s\n', ...
            'rank', 'real(lambda)', 'imag(lambda)', 'rel_res', 'abs_res', ...
            'matched_true', 'err_true');

        for t = 1:numel(idx_top)
            j = idx_top(t);
            [err_to_true, idx_true] = min(abs(lam_tls_all(j) - lam_true(:)));
            matched_true = lam_true(idx_true);
            fprintf('%4d  %16.12f  %11.3e  %12.3e  %12.3e  %16.12f  %12.3e\n', ...
                t, real(lam_tls_all(j)), imag(lam_tls_all(j)), rel_res(j), ...
                abs_res(j), matched_true, err_to_true);
        end
    end
end

function QB = build_QB_polygon(geom, cfg, lam)
% Centralized place for:
% 1. build A(lambda)
% 2. optional column normalization
% 3. QR truncation
% 4. return boundary block QB

    A = build_A_polygon(geom.P, geom.corner_list, geom.Mcorner, lam);

    if isfield(cfg, 'normalize_columns') && cfg.normalize_columns
        A = A ./ max(vecnorm(A), 1e-300);
    end

    [Q, R, ~] = qr(A, 0);

    if isfield(cfg, 'qr_tau') && ~isempty(cfg.qr_tau)
        d = abs(diag(R));
        if ~isempty(d)
            r = find(d >= cfg.qr_tau * d(1), 1, 'last');
            if isempty(r)
                Q = Q(:, []);
            else
                Q = Q(:, 1:r);
            end
        end
    end

    QB = Q(1:geom.mB, :);
end
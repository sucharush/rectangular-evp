function [A, B] = build_barycentric_star_pencil(D, z_nodes, w_weights)
% Star-topology linear pencil strictly conditioned for TLS eigenvalue extraction
% (direct barycentric-form linearization).

    d = numel(z_nodes) - 1;
    [p, q] = size(D{1});

    nrows = p + (d+1)*q;
    ncols = (d+2)*q;

    A = sparse(nrows, ncols);
    B = sparse(nrows, ncols);

    % --- CRITICAL FIX 1: Normalize Barycentric Weights ---
    % This eliminates internal dynamic range explosions without altering the roots.
    w_weights = w_weights / max(abs(w_weights));

    % 1. Top block row (Physics equation)
    for j = 1:(d+1)
        cols = j*q + (1:q);
        A(1:p, cols) = w_weights(j) * D{j};
    end

    % 2. Lower block rows (Math recurrence: x + z_i * y_i)
    Iq = speye(q);
    for i = 1:(d+1)
        rows      = p + (i-1)*q + (1:q);
        col_x     = 1:q;
        col_yi    = i*q + (1:q);

        A(rows, col_x)  = 1 * Iq;
        A(rows, col_yi) = z_nodes(i) * Iq;
        B(rows, col_yi) = 1 * Iq;
    end

    % --- CRITICAL FIX 2: Deliberate TLS Weighting ---
    norm_top = norm(A(1:p, :), 'inf');
    norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');

    ratio = norm_bot / max(norm_top, eps);
    if ratio > 1e2 || ratio < 1e-2
        gamma = sqrt(ratio);
        A(1:p, :) = gamma * A(1:p, :);
        B(1:p, :) = gamma * B(1:p, :);
    end
end

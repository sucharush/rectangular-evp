function [A, B] = build_newton_rect_pencil(D, sigma, beta, h, k)
% Direct rectangular analogue of Theorem 3 (Newton-form linearization).
%
% D{1},...,D{m+1} correspond to D_0,...,D_m, each p-by-q
% sigma, beta, h, k are length-m

    m = numel(beta);
    [p,q] = size(D{1});

    nrows = p + (m-1)*q;
    ncols = m*q;

    A = sparse(nrows, ncols);
    B = sparse(nrows, ncols);

    hm = h(m);
    km = k(m);
    betam = beta(m);

    % top block row
    for j = 1:(m-1)
        cols = (j-1)*q + (1:q);
        A(1:p, cols) = hm * D{j};
        B(1:p, cols) = km * D{j};
    end

    cols = (m-1)*q + (1:q);
    % FIX 1: Removed spurious 'hm' from the D{m+1} tail correction
    A(1:p, cols) = hm * D{m} - (sigma(m) / betam) * D{m+1};
    B(1:p, cols) = km * D{m} - (1 / betam) * D{m+1};

    Iq = speye(q);

    % lower block rows
    for i = 1:(m-1)
        rows      = p + (i-1)*q + (1:q);
        col_left  = (i-1)*q + (1:q);
        col_right = i*q     + (1:q);

        % no 'h(i)' from the last term
        A(rows, col_left)  = sigma(i) * Iq;
        A(rows, col_right) = h(i) * beta(i)  * Iq;

        B(rows, col_left)  = 1 * Iq;
        B(rows, col_right) = k(i) * beta(i)  * Iq;
    end
    % Compute a rough norm of the top row vs the identity blocks
    norm_top = norm(A(1:p, :), 'inf') + norm(B(1:p, :), 'inf');
    norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');

    % Scale the top row equations to match the lower recurrence equations
    gamma = norm_bot / max(norm_top, 1e-14);
    A(1:p, :) = gamma * A(1:p, :);
    B(1:p, :) = gamma * B(1:p, :);
end

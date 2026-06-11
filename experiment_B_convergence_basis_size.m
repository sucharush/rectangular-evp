clear; close all;
addpath(fileparts(mfilename('fullpath')));

alpha = 3;
beta = 0.04;
restol = 1e-1;

n_low = 15:2:40;
n_mid = 24:2:60;

lam_true_low = true_mode(alpha, beta, 7);
lam_true_mid = true_mode(alpha, beta, 12);

err_low_square = zeros(size(n_low));
err_low_rect = zeros(size(n_low));

for jj = 1:numel(n_low)
    n = n_low(jj);
    err_low_square(jj) = mode_error(n, n, alpha, beta, 0, 20, restol, lam_true_low);
    err_low_rect(jj) = mode_error(n, n+10, alpha, beta, 0, 20, restol, lam_true_low);
end

err_mid_square = zeros(size(n_mid));
err_mid_rect = zeros(size(n_mid));

for jj = 1:numel(n_mid)
    n = n_mid(jj);
    err_mid_square(jj) = mode_error(n, n, alpha, beta, 0, 25, restol, lam_true_mid);
    err_mid_rect(jj) = mode_error(n, ceil(1.5*n), alpha, beta, 0, 25, restol, lam_true_mid);
    err_mid_rect(jj) = mode_error(n, n+10, alpha, beta, 0, 25, restol, lam_true_mid);
end

figure;
hold on;
h1 = plot(n_low, err_low_square, 'k--', 'LineWidth', 1.2);
h2 = plot(n_low, err_low_rect, 'b-', 'LineWidth', 1.2);
hold off;
xlabel('n');
ylabel('|\lambda - \lambda^*|');
legend([h1, h2], 'm = n', 'm = ceil(2 n)', 'Location', 'best');
grid on;
box on;
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

figure;
hold on;
h1 = plot(n_mid, err_mid_square, 'k--', 'LineWidth', 1.2);
h2 = plot(n_mid, err_mid_rect, 'b-', 'LineWidth', 1.2);
hold off;
xlabel('n');
ylabel('|\lambda - \lambda^*|');
legend([h1, h2], 'm = n', 'm = ceil(1.5 n)', 'Location', 'best');
grid on;
box on;
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

function err = mode_error(n, m, alpha, beta, minbound, maxbound, restol, lam_true)
    coeffs = build_discretization(n, m, alpha, beta);
    [A, B, ~, ~, F_res] = poly_rect_linearization(coeffs);
    [vec_tls_all, lam_tls_all] = tls_pencil_eigs(A, B);

    keep = isfinite(real(lam_tls_all)) & isfinite(imag(lam_tls_all));
    lam_tls_all = lam_tls_all(keep);
    vec_tls_all = vec_tls_all(:, keep);

    if m == n
        keep_comp = candidate_window(lam_tls_all, minbound, maxbound);
        lam_comp = lam_tls_all(keep_comp);
    else
        [lam_comp, ~, ~, ~, ~] = tls_refine( ...
            lam_tls_all, vec_tls_all, F_res, lam_true, minbound, maxbound, restol);
    end

    if isempty(lam_comp)
        err = Inf;
        return;
    end

    err = min(abs(lam_comp - lam_true));
end

function lam_true = true_mode(alpha, beta, target_im)
    Kmax = max(200, ceil(2 * target_im / pi) + 8);
    k = (1:Kmax).';
    mu = (k * pi / 2).^2;
    damp_true = 2 * alpha + beta * mu;
    disc_true = damp_true.^2 - 4 * mu;
    lam_true_pos = (-damp_true + sqrt(disc_true)) / 2;
    lam_true_neg = (-damp_true - sqrt(disc_true)) / 2;
    lam_true_all = [lam_true_pos; lam_true_neg];
    lam_true_all = lam_true_all(imag(lam_true_all) >= 0);
    [~, idx] = min(abs(imag(lam_true_all) - target_im));
    lam_true = lam_true_all(idx);
end

function [coeffs, x, mats] = build_discretization(n, m, alpha, beta)
    if nargin < 3
        alpha = 1.0;
    end
    if nargin < 4
        beta = 0.01;
    end

    if m < n
        error('Require m >= n. Use m = n for square collocation and m > n for oversampling.');
    end

    ii = (1:m).';
    x = -cos(ii * pi / (m + 1));

    T = zeros(m, n);
    T(:, 1) = 1;

    if n >= 2
        T(:, 2) = x;
        for j = 2:n-1
            T(:, j+1) = 2 * x .* T(:, j) - T(:, j-1);
        end
    end

    Ujm1 = zeros(m, n);

    if n >= 2
        U0 = ones(m, 1);
        Ujm1(:, 2) = U0;

        if n >= 3
            U1 = 2 * x;
            Ujm1(:, 3) = U1;

            Uprev = U0;
            Ucurr = U1;
            for j = 3:n-1
                Unext = 2 * x .* Ucurr - Uprev;
                Ujm1(:, j+1) = Unext;
                Uprev = Ucurr;
                Ucurr = Unext;
            end
        end
    end

    M = zeros(m, n);
    K = zeros(m, n);

    for j = 0:n-1
        col = j + 1;
        M(:, col) = (1 - x.^2) .* T(:, col);
        K(:, col) = 3 * j * x .* Ujm1(:, col) + (j^2 + 2) * T(:, col);
    end

    C = 2 * alpha * M + beta * K;
    coeffs = {K, C, M};

    if nargout > 2
        mats.T = T;
        mats.Ujm1 = Ujm1;
        mats.M = M;
        mats.C = C;
        mats.K = K;
    end
end

function [A, B, F, F_mv, F_res] = poly_rect_linearization(Acell)
    N = numel(Acell) - 1;
    if N < 1
        error('Acell must contain at least {A0, A1}.');
    end

    [n, m] = size(Acell{1});

    for k = 2:numel(Acell)
        if ~isequal(size(Acell{k}), [n, m])
            error('All Ai must have the same size n-by-m.');
        end
    end

    rows = n + (N - 1) * m;
    cols = N * m;

    A = zeros(rows, cols);
    B = zeros(rows, cols);

    for j = 1:N
        col_idx = (j - 1) * m + (1:m);
        A(1:n, col_idx) = Acell{N - j + 1};
    end

    B(1:n, 1:m) = -Acell{N + 1};

    for i = 1:N-1
        row_idx = n + (i - 1) * m + (1:m);
        col_left = (i - 1) * m + (1:m);
        col_right = i * m + (1:m);

        A(row_idx, col_left) = eye(m);
        B(row_idx, col_right) = eye(m);
    end

    function mpoly = build_matrix_poly(lam)
        mpoly = zeros(n, m);
        for i = 1:N+1
            mpoly = mpoly + (lam^(i - 1)) * Acell{i};
        end
    end

    function mvp = build_matvec_poly(lam, y)
        mvp = zeros(n, 1);
        for i = 1:N+1
            mvp = mvp + (lam^(i - 1)) * Acell{i} * y;
        end
    end

    function [res_rel, res_abs, denom] = build_residual(lam, y)
        mvp = zeros(n, 1);
        denom = 0;

        for i = 1:N+1
            term = (lam^(i - 1)) * (Acell{i} * y);
            mvp = mvp + term;
            denom = denom + norm(term);
        end

        res_abs = norm(mvp);

        if denom == 0
            res_rel = Inf;
        else
            res_rel = res_abs / denom;
        end
    end

    F = @build_matrix_poly;
    F_mv = @build_matvec_poly;
    F_res = @build_residual;
end

function [lam_refined, vec_refined, res_refined, lam_exact, err_exact] = ...
    tls_refine(lam_tls_all, vec_tls_all, F_res, lam_true_all, minbound, maxbound, restol)
    n = size(vec_tls_all, 1) / 2;
    m = numel(lam_tls_all);
    lam_refined = [];
    vec_refined = zeros(n, 0);
    res_refined = [];
    lam_exact = [];
    err_exact = [];
    for ii = 1:m
        lam = lam_tls_all(ii);
        vec = vec_tls_all(:, ii);
        y1 = vec(n + 1:end);
        y2 = vec(1:n) / lam;
        if length(y1) ~= length(y2)
            error('extracted eigenvalues have different dimension');
        end
        U = [y1, y2];
        [Q, ~, ~] = svd(U, 'econ');
        y = Q(:, 1);
        res = F_res(lam, y);
        if ~candidate_window(lam, minbound, maxbound) || res > restol
            continue;
        end
        [err, jj] = min(abs(lam - lam_true_all));
        lam_refined(end+1, 1) = lam; %#ok<AGROW>
        vec_refined(:, end+1) = y; %#ok<AGROW>
        res_refined(end+1, 1) = res; %#ok<AGROW>
        lam_exact(end+1, 1) = lam_true_all(jj); %#ok<AGROW>
        err_exact(end+1, 1) = err; %#ok<AGROW>
    end
end

function keep = candidate_window(lam, minbound, maxbound)
    keep = imag(lam) >= 0 ...
        & abs(real(lam)) <= maxbound ...
        & abs(imag(lam)) <= maxbound ...
        & (abs(real(lam)) >= minbound | abs(imag(lam)) >= minbound);
end

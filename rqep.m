% clear; close all;

n = 15;
m_list = [n, ceil(n+10)];
alpha = 3;
beta = 0.04;
minbound = 1;
maxbound = 12;
restol = 1e-3;

Kmax = max(n + 8, ceil(2 * maxbound / pi) + 8);
k = (1:Kmax).';
mu = (k * pi / 2).^2;
damp_true = 2 * alpha + beta * mu;
disc_true = damp_true.^2 - 4 * mu;
lam_true_pos = (-damp_true + sqrt(disc_true)) / 2;
lam_true_neg = (-damp_true - sqrt(disc_true)) / 2;
lam_true_all = [lam_true_pos; lam_true_neg];
lam_true_plot = lam_true_all( ...
    abs(real(lam_true_all)) <= maxbound ...
    & abs(imag(lam_true_all)) <= maxbound);

for jj = 1:numel(m_list)
    m = m_list(jj);

    coeffs = build_discretization(n, m, alpha, beta);
    [A, B, ~, ~, F_res] = poly_rect_linearization(coeffs);
    [vec_tls_all, lam_tls_all] = tls_pencil_eigs(A, B);

    keep = isfinite(real(lam_tls_all)) & isfinite(imag(lam_tls_all));
    lam_tls_all = lam_tls_all(keep);
    vec_tls_all = vec_tls_all(:, keep);

    keep_plot = abs(real(lam_tls_all)) <= maxbound ...
        & abs(imag(lam_tls_all)) <= maxbound;
    lam_tls_plot = lam_tls_all(keep_plot);
    vec_tls_plot = vec_tls_all(:, keep_plot);

    figure;
    plot(real(lam_true_plot), imag(lam_true_plot), 'ko', 'MarkerSize', 6, 'LineWidth', 1.2);
    hold on;

    if m == n
        plot(real(lam_tls_plot), imag(lam_tls_plot), 'bx', 'MarkerSize', 7, 'LineWidth', 1.2);
        legend('exact', 'computed', 'Location', 'best');
    else
        [lam_accept, ~, res_accept, lam_exact, err_exact] = tls_refine( ...
            lam_tls_plot, vec_tls_plot, F_res, lam_true_plot, restol);
        keep_table = selected_table_window(lam_accept, minbound, maxbound);
        match_table = table( ...
            lam_accept(keep_table), lam_exact(keep_table), ...
            res_accept(keep_table), err_exact(keep_table), ...
            'VariableNames', {'lam_accept', 'lam_exact', 'res_accept', 'err_exact'});
        fprintf('m = %d accepted eigenvalue error table:\n', m);
        disp(match_table);

        keep_reject = true(size(lam_tls_plot));
        for ii = 1:numel(lam_accept)
            [~, idx] = min(abs(lam_tls_plot - lam_accept(ii)) + (~keep_reject) * 1e100);
            keep_reject(idx) = false;
        end
        lam_reject = lam_tls_plot(keep_reject);

        plot(real(lam_accept), imag(lam_accept), 'bx', 'MarkerSize', 7, 'LineWidth', 1.2);
        plot(real(lam_reject), imag(lam_reject), 'r+', 'MarkerSize', 7, 'LineWidth', 1.2);
        legend('exact', 'accepted', 'rejected', 'Location', 'best', 'FontSize', 13);
    end

    hold off;
    xlabel('Re(\lambda)');
    ylabel('Im(\lambda)');
    box on;
    xlim([-7, 0]);
    ylim([-maxbound, maxbound]);
    grid on;
    save_plot_eps("qep_vals")
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
            % res_rel = res_abs / denom;
            res_rel = res_abs;
        end
    end

    F = @build_matrix_poly;
    F_mv = @build_matvec_poly;
    F_res = @build_residual;
end

function [lam_refined, vec_refined, res_refined, lam_exact, err_exact] = ...
    tls_refine(lam_tls_all, vec_tls_all, F_res, lam_true_all, restol)
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
        if res > restol
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

function keep = selected_table_window(lam, minbound, maxbound)
    keep = imag(lam) >= 0 ...
        & abs(real(lam)) <= maxbound ...
        & abs(imag(lam)) <= maxbound ...
        & (abs(real(lam)) >= minbound | abs(imag(lam)) >= minbound);
end
function save_plot_eps(name)
    plot_dir = 'saved_plots';
    if exist(plot_dir, 'dir') ~= 7
        warning('save_plot_eps:MissingDir', ...
            'Directory "%s" does not exist. Skip saving plot "%s".', plot_dir, name);
        return;
    end

    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6.6, 5.8]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPositionMode', 'auto');

    filename = fullfile(plot_dir, sprintf('%s.eps', name));
    print(gcf, filename, '-depsc2', '-painters');
end

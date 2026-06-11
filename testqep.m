clear; close all;

n = 30;
m = 30;
alpha = 3;
beta = 0.04;
minbound = 5;
maxbound = 15;
restol = 1e-5;

coeffs = build_discretization(n, m, alpha, beta);
[A, B, ~, F_mv, F_res] = poly_rect_linearization(coeffs);

% ------------------------------------------------------------
% true eigenvalues for
%   lambda^2 + (2*alpha + beta*mu_k)*lambda + mu_k = 0
% ------------------------------------------------------------
Kmax = max(n + 8, ceil(2 * maxbound / pi) + 8);
k = (1:Kmax).';
mu = (k * pi / 2).^2;
damp_true = 2 * alpha + beta * mu;
disc_true = damp_true.^2 - 4 * mu;
lam_true_pos = (-damp_true + sqrt(disc_true)) / 2;
lam_true_neg = (-damp_true - sqrt(disc_true)) / 2;
lam_true_all = [lam_true_pos; lam_true_neg];

keep_true = imag(lam_true_all) >= 0 ...
    & abs(real(lam_true_all)) <= maxbound ...
    & abs(imag(lam_true_all)) <= maxbound ...
    & (abs(real(lam_true_all)) >= minbound | abs(imag(lam_true_all)) >= minbound);
lam_true_all = lam_true_all(keep_true);

% ------------------------------------------------------------
% TLS eigenvalues and error table
% ------------------------------------------------------------
[vec_tls_all, lam_tls_all] = tls_pencil_eigs(A, B);

keep = isfinite(real(lam_tls_all)) & isfinite(imag(lam_tls_all));
lam_tls_all = lam_tls_all(keep);
vec_tls_all = vec_tls_all(:, keep);

[lam_refined, ~, res_refined, lam_exact, err_exact] = tls_refine( ...
    lam_tls_all, vec_tls_all, F_mv, F_res, lam_true_all, minbound, maxbound, restol);

match_table = table(lam_refined, lam_exact, res_refined, err_exact);
disp(match_table);

function [coeffs, x, mats] = build_discretization(n, m, alpha, beta)
%BUILD_DISCRETIZATION
% Rectangular QEP for
%
%   lambda^2 v + 2*alpha*lambda v - beta*lambda*v'' - v'' = 0,
%
% with Dirichlet BC absorbed into the basis
%
%   phi_j(x) = (1 - x^2) T_j(x),   j = 0,...,n-1.
%
% Output:
%   coeffs = {K, C, M}
% so that
%   (K + lambda*C + lambda^2*M)c is approximately 0.
%
% Here
%   M_ij = phi_j(x_i),
%   C_ij = 2*alpha*phi_j(x_i) - beta*phi_j''(x_i),
%   K_ij = -phi_j''(x_i).

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
% Linearize
%   F(lambda) = A0 + lambda A1 + ... + lambda^N AN
% into
%   A y - lambda B y = 0
% with
%   y = [lambda^(N-1)x; lambda^(N-2)x; ...; x].

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
    tls_refine(lam_tls_all, vec_tls_all, F_mv, F_res, lam_true_all, minbound, maxbound, restol)
%TLS_REFINE Postprocess TLS candidates for the rectangular QEP.

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

        in_window = imag(lam) >= 0 ...
            && abs(real(lam)) <= maxbound ...
            && abs(imag(lam)) <= maxbound ...
            && (abs(real(lam)) >= minbound && abs(imag(lam)) >= minbound);
        if ~in_window
            continue;
        end

        top_block = vec(1:n);
        bottom_block = vec(n+1:end);

        if length(top_block) ~= length(bottom_block)
            error('Extracted eigenvector blocks have different dimensions.');
        end

        if abs(lam) < eps
            S = bottom_block;
        else
            S = [bottom_block, top_block / lam];
        end

        [U, R] = qr(S, 0);
        diagR = abs(diag(R));
        if isempty(diagR)
            continue;
        end

        keep = diagR > 1e-12 * max(diagR);
        U = U(:, keep);

        if isempty(U)
            continue;
        end

        FU = zeros(numel(F_mv(lam, U(:, 1))), size(U, 2));
        for jj = 1:size(U, 2)
            FU(:, jj) = F_mv(lam, U(:, jj));
        end

        [~, ~, Vsmall] = svd(FU, 'econ');
        z = Vsmall(:, end);
        y = U * z;
        y = y / norm(y);

        res = F_res(lam, y);

        if res > restol
            continue;
        end

        [err, jj_true] = min(abs(lam - lam_true_all));

        lam_refined(end+1, 1) = lam; %#ok<AGROW>
        vec_refined(:, end+1) = y; %#ok<AGROW>
        res_refined(end+1, 1) = res; %#ok<AGROW>
        lam_exact(end+1, 1) = lam_true_all(jj_true); %#ok<AGROW>
        err_exact(end+1, 1) = err; %#ok<AGROW>
    end
end

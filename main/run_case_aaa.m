clear; close all; clc;
rng(0);

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(fullfile(project_root, 'cases'));
addpath(fullfile(project_root, 'polygon'));
addpath(fullfile(project_root, 'problem_builders'));
addpath(fullfile(project_root, 'solvers'));
addpath(fullfile(project_root, 'rational'));
addpath(fullfile(project_root, 'core'));

% ============================================================
% Example 1: polygon backend
% ============================================================
cfg = case_polygon_drum('left');
cfg.qr_tau = 0;      % no truncation, keep QB size fixed
cfg.nb_per_edge = 10;
cfg.Mcorner = 10;
cfg.nI = 50;

problem = build_polygon_problem(cfg);

A_op  = problem.ops.A;
mB    = problem.meta.mB;
QB_op = @(lam) build_QB_from_Aop(A_op, mB, lam);

% ============================================================
% Example 2: synthetic operator backend
% Uncomment if needed
% ============================================================
% mB = 120;
% mI = 50;
% ncol = 60;
% [QB_op, A_op] = make_QB_op_poly(mB, mI, ncol, 1);

a = 7;
b = 11;

% ============================================================
% AAA options (base)
% ============================================================
opts = struct();
opts.interval = [a, b];

opts.aaa = struct();
opts.aaa.nCand = 1000;
opts.aaa.aaa_tol = 1e-11;
opts.aaa.mmax = 80;
opts.aaa.seed = 0;
opts.aaa.max_norm = true;
opts.aaa.svd_update = true;
opts.aaa.beta_mode = 'ones';

opts.verify = struct();
opts.verify.compute_bary_error = true;
%%
% ============================================================
% Run for ell = 2 and 4, average over k runs
% Keep the final stored result at seed = 0
% ============================================================
ell_list = [2, 4];
results_direct = cell(numel(ell_list), 1);
results_proc   = cell(numel(ell_list), 1);

k_runs = 1; %by default
seed_list = [1:(k_runs-1), 0];   % last run is seed = 0

avg_direct_degree = zeros(numel(ell_list), 1);
avg_proc_degree   = zeros(numel(ell_list), 1);
avg_direct_error  = zeros(numel(ell_list), 1);
avg_proc_error    = zeros(numel(ell_list), 1);
avg_direct_time = zeros(numel(ell_list), 1);
avg_proc_time = zeros(numel(ell_list), 1);

for ii = 1:numel(ell_list)
    ell = ell_list(ii);

    fprintf('\n================ ell = %d ================\n', ell);

    sum_direct_degree = 0;
    sum_proc_degree   = 0;
    sum_direct_error  = 0;
    sum_proc_error    = 0;
    sum_direct_time = 0;
    sum_proc_time = 0;

    for rr = 1:k_runs
        this_seed = seed_list(rr);

        fprintf('  run %d / %d, seed = %d\n', rr, k_runs, this_seed);

        % ----- direct mode -----
        opts_direct = opts;
        opts_direct.aaa.method = 'direct';
        opts_direct.aaa.ell = ell;
        opts_direct.aaa.seed = this_seed;

        result_direct_tmp = solve_aaa(QB_op, opts_direct);

        % ----- Procrustes-from-A mode -----
        opts_proc = opts;
        opts_proc.aaa.method = 'procrustes_from_a';
        opts_proc.aaa.ell = ell;
        opts_proc.aaa.seed = this_seed;
        opts_proc.aaa.mB = mB;
        opts_proc.aaa.sign_fix = true;
        opts_proc.aaa.store_Q  = true;

        result_proc_tmp = solve_aaa(A_op, opts_proc);

        % ----- sanity check -----
        Z_direct = result_direct_tmp.raw.Z(:);
        Z_proc   = result_proc_tmp.raw.Z(:);

        if numel(Z_direct) ~= numel(Z_proc) || any(abs(Z_direct - Z_proc) > 1e-12)
            error('The sampled Z grids are not identical for ell = %d, seed = %d.', ell, this_seed);
        end

        % ----- accumulate averages -----
        sum_direct_error  = sum_direct_error  + result_direct_tmp.verification.mean_relerr_bary;
        sum_proc_error    = sum_proc_error    + result_proc_tmp.verification.mean_relerr_bary;
        sum_direct_degree = sum_direct_degree + result_direct_tmp.raw.m;
        sum_proc_degree   = sum_proc_degree   + result_proc_tmp.raw.m;
        sum_direct_time = sum_direct_time + result_direct_tmp.raw.elapsed_time;
        sum_proc_time = sum_proc_time + result_proc_tmp.raw.elapsed_time;

        % ----- keep final run (seed = 0) under original names -----
        if rr == k_runs
            results_direct{ii} = result_direct_tmp;
            results_proc{ii}   = result_proc_tmp;
        end
    end

    avg_direct_error(ii)  = sum_direct_error  / k_runs;
    avg_proc_error(ii)    = sum_proc_error    / k_runs;
    avg_direct_degree(ii) = sum_direct_degree / k_runs;
    avg_proc_degree(ii)   = sum_proc_degree   / k_runs;
    avg_direct_time(ii) = sum_direct_time / k_runs;
    avg_proc_time(ii) = sum_proc_time / k_runs;

    fprintf('Average barycentric rel. Frobenius error (normalized)     : %.6e\n', ...
        avg_direct_error(ii));
    fprintf('Average barycentric rel. Frobenius error (Procrustes) : %.6e\n', ...
        avg_proc_error(ii));
    fprintf('Average AAA degree / support points (normalized)     : %.6f\n', ...
        avg_direct_degree(ii));
    fprintf('Average AAA degree / support points (Procrustes) : %.6f\n', ...
        avg_proc_degree(ii));
    fprintf('Average AAA time (normalized) : %.6f\n', ...
        avg_direct_time(ii));
    fprintf('Average AAA time (Procrustes) : %.6f\n', ...
        avg_proc_time(ii));

    fprintf('Stored in results_* for plotting: final run with seed = 0\n');
end
%%
% ============================================================
% Plot 1: stacked view for ell = 2 and 4
% ============================================================
figure;

t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

styles = struct();
styles.direct = '-';
styles.proc   = '--';

for ii = 1:numel(ell_list)
    ell = ell_list(ii);

    nexttile;
    ax = gca;
    hold(ax, 'on');
    set(ax, 'YScale', 'log');
    box(ax, 'on');

    Z_direct = results_direct{ii}.raw.Z(:);
    Z_proc   = results_proc{ii}.raw.Z(:);

    plot(ax, Z_direct, results_direct{ii}.verification.relerr_bary, ...
        styles.direct, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('raw, l=%d', ell));

    plot(ax, Z_proc, results_proc{ii}.verification.relerr_bary, ...
        styles.proc, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Procrustes, l=%d', ell));

    xlim([a, b]);
    grid on;
    ylabel('rel. Fro error');
    legend();
    % legend('Location', 'best');

    if ell == 2
        ylim([1e-10, 1e-7]);
        title('l = 2');
        set(gca, 'XTickLabel', []);
    elseif ell == 4
        ylim([1e-12, 1e-8]);
        title('l = 4');
        xlabel('z');
    end
end
if ~exist('saved_plots', 'dir')
    mkdir('saved_plots');
end
% title(t, 'Barycentric approximation error');
% export_fig("saved_plots/error_bary.eps")
saveas(gcf, 'saved_plots/error_bary.eps', 'epsc');
%%
% ============================================================
% Plot 2: consecutive change on [a,b]
% ============================================================
lam_grid = linspace(a, b, 1000);
nlam = numel(lam_grid);

A0 = A_op(lam_grid(1));
[Q0, R0] = qr(A0, 0);

s0 = sign(diag(R0));
s0(s0 == 0) = 1;
Q0 = Q0 * diag(s0);

r = size(Q0, 2);

QB_raw_all  = zeros(mB, r, nlam);
QB_proc_all = zeros(mB, r, nlam);
QB_proj_all = zeros(mB, mB, nlam);

QB_raw_all(:,:,1)  = Q0(1:mB, :);
QB_proc_all(:,:,1) = Q0(1:mB, :);
QB_proj_all(:,:,1) = Q0(1:mB, :) * Q0(1:mB, :)';

Q_prev = Q0;

for k = 2:nlam
    Ak = A_op(lam_grid(k));
    [Qk, Rk] = qr(Ak, 0);

    s = sign(diag(Rk));
    s(s == 0) = 1;
    Qk = Qk * diag(s);

    % orthogonal Procrustes: align Qk to Q_prev
    [U, ~, V] = svd(Q_prev' * Qk, 'econ');
    Ralign = V * U';
    Qk_proc = Qk * Ralign;

    QB_raw_all(:,:,k)  = Qk(1:mB, :);
    QB_proc_all(:,:,k) = Qk_proc(1:mB, :);
    QB_proj_all(:,:,k) = Qk(1:mB, :) * Qk(1:mB, :)';

    Q_prev = Qk_proc;
end

diff_raw  = zeros(nlam-1, 1);
diff_proc = zeros(nlam-1, 1);
diff_proj = zeros(nlam-1, 1);

for k = 2:nlam
    diff_raw(k-1)  = norm(QB_raw_all(:,:,k)  - QB_raw_all(:,:,k-1),  'fro');
    diff_proc(k-1) = norm(QB_proc_all(:,:,k) - QB_proc_all(:,:,k-1), 'fro');
    diff_proj(k-1) = norm(QB_proj_all(:,:,k) - QB_proj_all(:,:,k-1), 'fro');
end
%%
figure;
semilogy(lam_grid(2:end), diff_raw,  '-',  'LineWidth', 1.2); hold on;
semilogy(lam_grid(2:end), diff_proc, '-.', 'LineWidth', 1.2); 
semilogy(lam_grid(2:end), diff_proj, '--', 'LineWidth', 1.2);
ylim([1e-3, 3e-2]);
set(gca, 'YScale', 'log');
ax = gca;
disp(ax.YScale)
grid on;
xlabel('\lambda');
% ylabel('||F_k - F_{k-1}||_F');
% legend('Q_B(\lambda)', 'Q_B(\lambda) + Procrustes', 'Q_B(\lambda)Q_B(\lambda)^*', ...
%     'Location', 'best');
legend('raw', 'Procrustes', 'projector');
title('Consecutive change');
% export_fig("saved_plots/continue.eps")
saveas(gcf, 'saved_plots/continue.eps', 'epsc');

%%
% ============================================================
% Plot 3: sigma_min on linearized pencils for all 4 outputs
% ============================================================
% lam_true = [ ...
%     7.248077862494475, ...
%     9.209294998335231, ...
%     10.596985691456322];
lam_true = [ ...
    7.247948913733, ...
    9.208978724772, ...
    10.594985597455];
delta = 0.03;

lam_coarse = linspace(a, b, 41);
lam_near = [ ...
    lam_true(1)-delta, lam_true(1), lam_true(1)+delta, ...
    lam_true(2)-delta, lam_true(2), lam_true(2)+delta, ...
    lam_true(3)-delta, lam_true(3), lam_true(3)+delta ...
];

lam_test = unique([lam_coarse, lam_near]);
lam_test = sort(lam_test);
nt = numel(lam_test);

idx4 = find(ell_list == 4, 1);

if isempty(idx4)
    error('ell = 4 not found in ell_list.');
end

tls_imag_tol = 1e-1;
tls_top_k = 10;
tls_reports = cell(2, 1);

figure;
ax = axes;
hold(ax, 'on');
set(ax, 'YScale', 'log');
box(ax, 'on');

outs = {results_direct{idx4}.raw, results_proc{idx4}.raw};
names = {'normalized, l=4', 'Procrustes, l=4'};
styles = {'o-', 's--'};

for k = 1:2
    out = outs{k};
    [A, B] = build_rectangular_linear_pencil(out.D, out.sigma, out.beta, out.h, out.k);

    sigvals = zeros(nt,1);
    for j = 1:nt
        s = svd(full(A - lam_test(j) * B), 'econ');
        sigvals(j) = s(end);
    end

    semilogy(ax, lam_test, sigvals, styles{k}, 'LineWidth', 1.2, ...
        'DisplayName', names{k});

    for j = 1:numel(lam_true)
        s_true = svd(full(A - lam_true(j) * B), 'econ');
        semilogy(ax, lam_true(j), s_true(end), 'ro', ...
            'MarkerSize', 6, 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end

    tls_reports{k} = report_tls_pencil_pairs(A, B, lam_true, [a, b], ...
        tls_imag_tol, tls_top_k, names{k});
end

for j = 1:numel(lam_true)
    xline(ax, lam_true(j), '--k', 'HandleVisibility', 'off');
end

xlabel('\lambda');
ylabel('\sigma_{min}(M-\lambda B)');
xlim([a, b]);
% legend('Location', 'best');
legend();
title('Smallest singular value of the linearized pencils, l=4');
grid on;
saveas(gcf, 'saved_plots/aaa_linearized.eps', 'epsc');

%%
% ============================================================
% Plot 4: TLS block-based post-filter on the linearized pencils
% ============================================================
tls_block_top_k = 20;
tls_block_reports = cell(2, 1);

for k = 2:2
    out = outs{k};
    tls_block_reports{k} = report_tls_block_postfilter(out, lam_true, [a, b], ...
        tls_imag_tol, tls_block_top_k, names{k});
end

%%
% %%%%%%%%%%%%%%%%%%%%
%     helpers
% %%%%%%%%%%%%%%%%%%%%
function QB = build_QB_from_Aop(A_op, mB, lam)
    A = A_op(lam);
    [Q, R] = qr(A, 0);

    s = sign(diag(R));
    s(s == 0) = 1;
    Q = Q * diag(s);

    QB = Q(1:mB, :);
end
function [A, B] = build_rectangular_linear_pencil(D, sigma, beta, h, k)
% Direct rectangular analogue of Theorem 3.
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

function report = report_tls_pencil_pairs(A, B, lam_true, interval, imag_tol, top_k, label)
    [vec_tls_all, lam_tls_all] = tls_pencil_eigs(full(A), full(B));

    lam_tls_all = lam_tls_all(:);
    nlam = numel(lam_tls_all);

    if size(vec_tls_all, 2) ~= nlam
        error('report_tls_pencil_pairs: inconsistent TLS eigenvector/eigenvalue sizes.');
    end

    abs_res = nan(nlam, 1);
    rel_res = nan(nlam, 1);
    matched_true = nan(nlam, 1);
    err_to_true = nan(nlam, 1);

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

        [err_to_true(j), idx_true] = min(abs(lam - lam_true(:)));
        matched_true(j) = lam_true(idx_true);
    end

    keep_mask = isfinite(real(lam_tls_all)) ...
        & isfinite(imag(lam_tls_all)) ...
        & isfinite(rel_res) ...
        & real(lam_tls_all) >= interval(1) ...
        & real(lam_tls_all) <= interval(2) ...
        & abs(imag(lam_tls_all)) <= imag_tol;

    idx_keep = find(keep_mask);
    [~, order] = sort(rel_res(idx_keep), 'ascend');
    idx_ranked = idx_keep(order);
    idx_top = idx_ranked(1:min(top_k, numel(idx_ranked)));

    fprintf('\nTLS pencil candidates: %s\n', label);
    fprintf('  prefilter: real(lambda) in [%.6f, %.6f], |imag(lambda)| <= %.3g\n', ...
        interval(1), interval(2), imag_tol);
    fprintf('  kept %d of %d TLS pairs; reporting top %d by relative residual.\n', ...
        numel(idx_keep), nlam, numel(idx_top));

    if isempty(idx_top)
        fprintf('  no TLS pairs passed the prefilter.\n');
    else
        fprintf('%4s  %16s  %11s  %12s  %12s  %16s  %12s\n', ...
            'rank', 'real(lambda)', 'imag(lambda)', 'rel_res', 'abs_res', ...
            'matched_true', 'err_true');

        for t = 1:numel(idx_top)
            j = idx_top(t);
            fprintf('%4d  %16.12f  %11.3e  %12.3e  %12.3e  %16.12f  %12.3e\n', ...
                t, real(lam_tls_all(j)), imag(lam_tls_all(j)), rel_res(j), ...
                abs_res(j), matched_true(j), err_to_true(j));
        end
    end

    report = struct();
    report.label = label;
    report.interval = interval;
    report.imag_tol = imag_tol;
    report.n_total = nlam;
    report.n_kept = numel(idx_keep);
    report.keep_mask = keep_mask;
    report.lam_all = lam_tls_all;
    report.abs_res_all = abs_res;
    report.rel_res_all = rel_res;
    report.matched_true_all = matched_true;
    report.err_to_true_all = err_to_true;
    report.top_indices = idx_top;
    report.top_lam = lam_tls_all(idx_top);
    report.top_abs_res = abs_res(idx_top);
    report.top_rel_res = rel_res(idx_top);
    report.top_matched_true = matched_true(idx_top);
    report.top_err_to_true = err_to_true(idx_top);
end

function report = report_tls_block_postfilter(out, lam_true, interval, imag_tol, top_k, label)
    [vec_tls_all, lam_tls_all] = tls_pencil_eigs(full(out.A), full(out.B));

    lam_tls_all = lam_tls_all(:);
    nlam = numel(lam_tls_all);
    m = out.m;
    q = out.q;

    if size(vec_tls_all, 2) ~= nlam
        error('report_tls_block_postfilter: inconsistent TLS eigenvector/eigenvalue sizes.');
    end
    if size(vec_tls_all, 1) ~= m * q
        error('report_tls_block_postfilter: TLS eigenvector length does not match m*q.');
    end

    rel_pencil_res = nan(nlam, 1);
    rb_score = nan(nlam, 1);
    matched_true = nan(nlam, 1);
    err_to_true = nan(nlam, 1);
    stable_block = ones(nlam, 1);
    stable_scale = nan(nlam, 1);
    y_norm = nan(nlam, 1);

    normA = norm(out.A, 'fro');
    normB = norm(out.B, 'fro');

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

        pencil_res = norm((out.A - lam * out.B) * x);
        rel_pencil_res(j) = pencil_res / max((normA + abs(lam) * normB) * xnorm, 1e-14);

        [y, block_id, scale_val] = local_extract_first_block_y(x, q);

        if isempty(y)
            continue;
        end

        Rb = aaa_eval_matrix_barycentric(lam, out.zj, out.wj, out.D);
        % QB_eval = QB_op(lam);
        % Rb = QB_eval* QB_eval;
        y_norm(j) = norm(y);
        rb_score(j) = norm(Rb * y) / max(y_norm(j), 1e-14);
        stable_block(j) = block_id;
        stable_scale(j) = scale_val;

        [err_to_true(j), idx_true] = min(abs(lam - lam_true(:)));
        matched_true(j) = lam_true(idx_true);
    end

    keep_mask = isfinite(real(lam_tls_all)) ...
        & isfinite(imag(lam_tls_all)) ...
        & isfinite(rel_pencil_res) ...
        & isfinite(rb_score) ...
        & real(lam_tls_all) >= interval(1) ...
        & real(lam_tls_all) <= interval(2) ...
        & abs(imag(lam_tls_all)) <= imag_tol;

    idx_keep = find(keep_mask);
    [~, order] = sort(rb_score(idx_keep), 'ascend');
    idx_ranked = idx_keep(order);
    idx_top = idx_ranked(1:min(top_k, numel(idx_ranked)));

    fprintf('\nTLS block post-filter: %s\n', label);
    fprintf('  prefilter: real(lambda) in [%.6f, %.6f], |imag(lambda)| <= %.3g\n', ...
        interval(1), interval(2), imag_tol);
    fprintf('  kept %d of %d TLS pairs; reporting top %d by post-filter score.\n', ...
        numel(idx_keep), nlam, numel(idx_top));

    if isempty(idx_top)
        fprintf('  no TLS pairs passed the prefilter.\n');
    else
        fprintf('%4s  %16s  %11s  %12s  %12s  %6s  %12s  %16s  %12s\n', ...
            'rank', 'real(lambda)', 'imag(lambda)', 'Rb_score', 'pencil_res', ...
            'block', 'scale', 'matched_true', 'err_true');

        for t = 1:numel(idx_top)
            j = idx_top(t);
            fprintf('%4d  %16.12f  %11.3e  %12.3e  %12.3e  %6d  %12.3e  %16.12f  %12.3e\n', ...
                t, real(lam_tls_all(j)), imag(lam_tls_all(j)), rb_score(j), ...
                rel_pencil_res(j), stable_block(j), stable_scale(j), matched_true(j), ...
                err_to_true(j));
        end
    end

    report = struct();
    report.label = label;
    report.interval = interval;
    report.imag_tol = imag_tol;
    report.n_total = nlam;
    report.n_kept = numel(idx_keep);
    report.keep_mask = keep_mask;
    report.lam_all = lam_tls_all;
    report.vec_all = vec_tls_all;
    report.rel_pencil_res_all = rel_pencil_res;
    report.rb_score_all = rb_score;
    report.matched_true_all = matched_true;
    report.err_to_true_all = err_to_true;
    report.stable_block_all = stable_block;
    report.stable_scale_all = stable_scale;
    report.y_norm_all = y_norm;
    report.top_indices = idx_top;
    report.top_lam = lam_tls_all(idx_top);
    report.top_rb_score = rb_score(idx_top);
    report.top_rel_pencil_res = rel_pencil_res(idx_top);
    report.top_matched_true = matched_true(idx_top);
    report.top_err_to_true = err_to_true(idx_top);
    report.top_stable_block = stable_block(idx_top);
    report.top_stable_scale = stable_scale(idx_top);
end

function [y, block_id, scale_val] = local_extract_first_block_y(x, q)
    if numel(x) < q
        y = [];
        block_id = [];
        scale_val = [];
        return;
    end

    y = x(1:q);
    block_id = 1;
    scale_val = norm(y);

    if ~isfinite(scale_val) || scale_val == 0
        y = [];
        block_id = [];
        scale_val = [];
    end
end

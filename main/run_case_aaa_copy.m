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
cfg.nI = 10;

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
ell_list = [4];
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
lam_true = [ ...
    7.248077955423, ...
    9.209295403292, ...
    10.596986220784];
% 7.248077955423  9.209295403292 10.596986220784
% lam_true = [ ...
%     7.247948913733, ...
%     9.208978724772, ...
%     10.594985597455];
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

tls_imag_tol = 1e-2;
tls_top_k = 10;
pencil_data = cell(2, 1);

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
    [A, B] = build_star_barycentric_pencil(out.D, out.zj, out.wj);
    tls_data = collect_tls_pencil_data(A, B, [a, b], tls_imag_tol);
    pencil_data{k} = struct('A', A, 'B', B, 'tls', tls_data);

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

    report_tls_pencil_pairs(tls_data, lam_true, tls_top_k, names{k});
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

for k = 2:2
    out = outs{k};
    report_tls_block_postfilter(out, pencil_data{k}, lam_true, tls_block_top_k, names{k});
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
% function [A, B] = build_star_barycentric_pencil(D, z_nodes, w_weights)
% % Star-topology rectangular linearization for the barycentric form.
% %
% % State blocks are [x, y_1, ..., y_{d+1}] with x the global q-vector.
% 
%     d = numel(z_nodes) - 1;
%     [p, q] = size(D{1});
% 
%     nrows = p + (d + 1) * q;
%     ncols = (d + 2) * q;
% 
%     A = sparse(nrows, ncols);
%     B = sparse(nrows, ncols);
% 
%     % Top block row: weighted barycentric sum, shifted by one block for x.
%     for j = 1:(d + 1)
%         cols = j * q + (1:q);
%         A(1:p, cols) = w_weights(j) * D{j};
%     end
% 
%     % Lower block rows: x + z_j y_j = lambda y_j.
%     Iq = speye(q);
%     for i = 1:(d + 1)
%         rows = p + (i - 1) * q + (1:q);
%         col_x = 1:q;
%         col_yi = i * q + (1:q);
% 
%         A(rows, col_x) = Iq;
%         A(rows, col_yi) = z_nodes(i) * Iq;
% 
%         B(rows, col_yi) = Iq;
%     end
% 
%     norm_top = norm(A(1:p, :), 'inf');
%     norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');
% 
%     gamma = norm_bot / max(norm_top, 1e-14);
%     A(1:p, :) = gamma * A(1:p, :);
% end
function [A, B] = build_star_barycentric_pencil(D, z_nodes, w_weights)
% Star-topology linear pencil strictly conditioned for TLS eigenvalue extraction.

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
    % Instead of blind block balancing, we explicitly measure the norms.
    % Assuming D{j} are properly scaled, norm_top should already be reasonable now.
    norm_top = norm(A(1:p, :), 'inf'); 
    norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');
    
    % Only apply a gentle gamma if the scale disparity is still
    % disastrously large...
    ratio = norm_bot / max(norm_top, eps);
    if ratio > 1e2 || ratio < 1e-2
        % Only act if there is a severe scale mismatch
        gamma = sqrt(ratio); % A milder preconditioner to prevent dominating the SVD
        A(1:p, :) = gamma * A(1:p, :);
        B(1:p, :) = gamma * B(1:p, :); 
    end
end

% Legacy chain-style barycentric linearization kept for reference.
% function [A, B] = build_barycentric_linear_pencil(D, z_nodes, w_weights)
% % Rectangular linearization built directly from the barycentric form.
% %
% % D{1},...,D{d+1} correspond to sampled blocks F(z_j), each p-by-q.
% % z_nodes and w_weights are the barycentric support points and weights.
%
%     d = numel(z_nodes) - 1;
%     [p, q] = size(D{1});
%
%     nrows = p + d * q;
%     ncols = (d + 1) * q;
%
%     A = sparse(nrows, ncols);
%     B = sparse(nrows, ncols);
%
%     for j = 1:(d + 1)
%         cols = (j - 1) * q + (1:q);
%         A(1:p, cols) = w_weights(j) * D{j};
%     end
%
%     Iq = speye(q);
%     for i = 1:d
%         rows = p + (i - 1) * q + (1:q);
%         col_left = (i - 1) * q + (1:q);
%         col_right = i * q + (1:q);
%
%         A(rows, col_left) = z_nodes(i) * Iq;
%         A(rows, col_right) = -z_nodes(i + 1) * Iq;
%
%         B(rows, col_left) = Iq;
%         B(rows, col_right) = -Iq;
%     end
%
%     norm_top = norm(A(1:p, :), 'inf');
%     norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');
%
%     gamma = norm_bot / max(norm_top, 1e-14);
%     A(1:p, :) = gamma * A(1:p, :);
% end

% Legacy Newton-style linearization kept for reference.
% function [A, B] = build_rectangular_linear_pencil(D, sigma, beta, h, k)
% % Direct rectangular analogue of Theorem 3.
% %
% % D{1},...,D{m+1} correspond to D_0,...,D_m, each p-by-q
% % sigma, beta, h, k are length-m
%
%     m = numel(beta);
%     [p,q] = size(D{1});
%
%     nrows = p + (m-1)*q;
%     ncols = m*q;
%
%     A = sparse(nrows, ncols);
%     B = sparse(nrows, ncols);
%
%     hm = h(m);
%     km = k(m);
%     betam = beta(m);
%
%     for j = 1:(m-1)
%         cols = (j-1)*q + (1:q);
%         A(1:p, cols) = hm * D{j};
%         B(1:p, cols) = km * D{j};
%     end
%
%     cols = (m-1)*q + (1:q);
%     A(1:p, cols) = hm * D{m} - (sigma(m) / betam) * D{m+1};
%     B(1:p, cols) = km * D{m} - (1 / betam) * D{m+1};
%
%     Iq = speye(q);
%     for i = 1:(m-1)
%         rows      = p + (i-1)*q + (1:q);
%         col_left  = (i-1)*q + (1:q);
%         col_right = i*q     + (1:q);
%
%         A(rows, col_left)  = sigma(i) * Iq;
%         A(rows, col_right) = h(i) * beta(i)  * Iq;
%
%         B(rows, col_left)  = Iq;
%         B(rows, col_right) = k(i) * beta(i)  * Iq;
%     end
%
%     norm_top = norm(A(1:p, :), 'inf') + norm(B(1:p, :), 'inf');
%     norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');
%
%     gamma = norm_bot / max(norm_top, 1e-14);
%     A(1:p, :) = gamma * A(1:p, :);
%     B(1:p, :) = gamma * B(1:p, :);
% end

function tls_data = collect_tls_pencil_data(A, B, interval, imag_tol)
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

function report_tls_pencil_pairs(tls_data, lam_true, top_k, label)
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

function report_tls_block_postfilter(out, pencil_data, lam_true, top_k, label)
    A = pencil_data.A;
    B = pencil_data.B;
    tls_data = pencil_data.tls;
    lam_tls_all = tls_data.lam_all;
    vec_tls_all = tls_data.vec_all;
    idx_keep = tls_data.idx_keep;

    nblocks = numel(out.zj);
    q = out.q;

    if size(vec_tls_all, 1) ~= (nblocks + 1) * q
        error('report_tls_block_postfilter: TLS eigenvector length does not match (numel(zj)+1)*q.');
    end

    nkeep = numel(idx_keep);
    reduced_score = nan(nkeep, 1);
    first_block_score = nan(nkeep, 1);
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

        y_first = local_extract_first_star_block(x, q);
        if ~isempty(y_first)
            y_first_norm = norm(y_first);
            first_block_score(t) = norm(Rb * y_first) / max(y_first_norm, 1e-14);
        end

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
        fprintf('%4s  %16s  %11s  %12s  %12s  %12s  %6s  %6s  %16s  %12s\n', ...
            'rank', 'real(lambda)', 'imag(lambda)', 'red_score', 'first_score', ...
            'pencil_res', 'dimU', 'nblk', 'matched_true', 'err_true');

        for t = 1:numel(idx_top)
            j = idx_top(t);
            fprintf('%4d  %16.12f  %11.3e  %12.3e  %12.3e  %12.3e  %6d  %6d  %16.12f  %12.3e\n', ...
                t, real(lam_keep(j)), imag(lam_keep(j)), reduced_score(j), ...
                first_block_score(j), rel_pencil_res(j), round(subspace_dim(j)), ...
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

    % Each y_j block provides an additional x candidate via x = (lambda - z_j) y_j.
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

function y = local_extract_first_star_block(x, q)
    if numel(x) < q
        y = [];
        return;
    end

    y = x(1:q);
    y_norm = norm(y);
    if ~isfinite(y_norm) || y_norm == 0
        y = [];
    end
end

clear; close all; clc;
seed = 0;
rng(seed);
last_seed = seed;

% ------------------------------------------------------------
% Linearization of the AAA rational surrogate:
%   'newton'      -> Newton-form rectangular pencil  (run_case_aaa_newton)
%   'barycentric' -> star-barycentric pencil         (run_case_aaa_bary)
% ------------------------------------------------------------
linearization = 'barycentric';

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(fullfile(project_root, 'cases'));
addpath(fullfile(project_root, 'polygon'));
addpath(fullfile(project_root, 'solvers'));
addpath(fullfile(project_root, 'rational'));
addpath(fullfile(project_root, 'core'));
addpath(fullfile(project_root, 'plots'));

strat = get_linearization_strategy(linearization);

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

k_runs = 5; %by default
seed_list = [1:(k_runs-1), last_seed];   % last run is seed = 0

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
        'DisplayName', sprintf('normalized, l=%d', ell));

    plot(ax, Z_proc, results_proc{ii}.verification.relerr_bary, ...
        styles.proc, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Procrustes, l=%d', ell));

    xlim([a, b]);
    grid on;
    legend(ax, 'FontSize', 13);

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
legend('normalized', 'Procrustes', 'projector', 'FontSize', 13);
saveas(gcf, 'saved_plots/continue.eps', 'epsc');

%%
% ============================================================
% Plot 3: sigma_min on linearized pencils for the two ell=4 outputs
% ============================================================
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

tls_imag_tol = strat.imag_tol;
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
    [A, B] = strat.build_pencil(out);
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
legend('FontSize', 13);
grid on;
saveas(gcf, sprintf('saved_plots/%s.eps', strat.save_name), 'epsc');

%%
% ============================================================
% Plot 4: TLS block-based post-filter on the linearized pencils
% ============================================================
tls_block_top_k = 20;

for k = 1:2
    out = outs{k};
    strat.postfilter(out, pencil_data{k}, lam_true, tls_block_top_k, names{k});
end

%%
% %%%%%%%%%%%%%%%%%%%%
%     helpers
% %%%%%%%%%%%%%%%%%%%%
function strat = get_linearization_strategy(name)
% Select the linearization-specific pieces: pencil builder, TLS imaginary
% tolerance, block post-filter, and the plot file name.
    switch lower(name)
        case 'newton'
            strat.name = 'newton';
            strat.build_pencil = @(out) build_newton_rect_pencil( ...
                out.D, out.sigma, out.beta, out.h, out.k);
            strat.imag_tol = 1e-1;
            strat.postfilter = @report_tls_newton_postfilter;
            strat.save_name = 'aaa_linearized_newton';
        case 'barycentric'
            strat.name = 'barycentric';
            strat.build_pencil = @(out) build_barycentric_star_pencil( ...
                out.D, out.zj, out.wj);
            strat.imag_tol = 1e-2;
            strat.postfilter = @report_tls_barycentric_postfilter;
            strat.save_name = 'aaa_linearized_bary';
        otherwise
            error('get_linearization_strategy: unknown linearization ''%s''.', name);
    end
end

function QB = build_QB_from_Aop(A_op, mB, lam)
    A = A_op(lam);
    [Q, R] = qr(A, 0);

    s = sign(diag(R));
    s(s == 0) = 1;
    Q = Q * diag(s);

    QB = Q(1:mB, :);
end

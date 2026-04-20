clear; close all; clc;
rng(1);

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(fullfile(project_root, 'cases'));
addpath(fullfile(project_root, 'polygon'));
addpath(fullfile(project_root, 'problem_builders'));
addpath(fullfile(project_root, 'solvers'));
addpath(fullfile(project_root, 'rational'));

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
opts.aaa.aaa_tol = 1e-10;
opts.aaa.mmax = 80;
opts.aaa.seed = 1;
opts.aaa.max_norm = true;
opts.aaa.svd_update = true;
opts.aaa.beta_mode = 'ones';

opts.verify = struct();
opts.verify.compute_bary_error = true;

% ============================================================
% Run for ell = 2 and 4
% ============================================================
ell_list = [2, 4];
results_direct = cell(numel(ell_list), 1);
results_proc   = cell(numel(ell_list), 1);

for ii = 1:numel(ell_list)
    ell = ell_list(ii);

    fprintf('\n================ ell = %d ================\n', ell);

    % ----- direct mode -----
    opts_direct = opts;
    opts_direct.aaa.method = 'direct';
    opts_direct.aaa.ell = ell;

    results_direct{ii} = solve_aaa(QB_op, opts_direct);

    % ----- Procrustes-from-A mode -----
    opts_proc = opts;
    opts_proc.aaa.method = 'procrustes_from_a';
    opts_proc.aaa.ell = ell;
    opts_proc.aaa.mB = mB;
    opts_proc.aaa.sign_fix = true;
    opts_proc.aaa.store_Q  = true;

    results_proc{ii} = solve_aaa(A_op, opts_proc);

    % ----- sanity check: same Z -----
    Z_direct = results_direct{ii}.raw.Z(:);
    Z_proc   = results_proc{ii}.raw.Z(:);

    if numel(Z_direct) ~= numel(Z_proc) || any(abs(Z_direct - Z_proc) > 1e-12)
        error('The sampled Z grids are not identical for ell = %d.', ell);
    end

    fprintf('Average barycentric rel. Frobenius error (direct)     : %.6e\n', ...
        results_direct{ii}.verification.mean_relerr_bary);
    fprintf('Average barycentric rel. Frobenius error (Procrustes) : %.6e\n', ...
        results_proc{ii}.verification.mean_relerr_bary);
    fprintf('AAA degree / support points (direct)     : %d\n', results_direct{ii}.raw.m);
    fprintf('AAA degree / support points (Procrustes) : %d\n', results_proc{ii}.raw.m);
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

    Z_direct = results_direct{ii}.raw.Z(:);
    Z_proc   = results_proc{ii}.raw.Z(:);

    plot(ax, Z_direct, results_direct{ii}.verification.relerr_bary, ...
        styles.direct, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('direct, l=%d', ell));

    plot(ax, Z_proc, results_proc{ii}.verification.relerr_bary, ...
        styles.proc, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Procrustes, l=%d', ell));

    xlim([a, b]);
    grid on;
    ylabel('rel. Fro error');
    legend('Location', 'best');

    if ell == 2
        ylim([1e-10, 1e-6]);
        title('l = 2');
        set(gca, 'XTickLabel', []);
    elseif ell == 4
        ylim([1e-11, 1e-7]);
        title('l = 4');
        xlabel('z');
    end
end

title(t, 'Barycentric approximation error');

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
ylabel('||F_k - F_{k-1}||_F');
legend('Q_B(\lambda)', 'Q_B(\lambda) + Procrustes', 'Q_B(\lambda)Q_B(\lambda)^*', ...
    'Location', 'best');
title('Consecutive change');

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
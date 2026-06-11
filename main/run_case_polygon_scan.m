clear; close all; clc;
rng(0);
% %% test polygon
% cfg = case_polygon_hshape();
% problem = build_polygon_problem(cfg);
% problem.ops.sigma(12.18)
% 
% %% test minimizer
% 
% obj = make_counted_objective(@(x) (x-2)^2 + 1);
% [lam_star, sig_star, info] = minimizer_hybrid(1, 3, obj, struct())
% 
% %% test scan
% cfg = case_polygon_hshape();
% problem = build_polygon_problem(cfg);
% 
% scan_opts = struct();
% scan_opts.lamvec = 10:0.01:13;
% scan_opts.detect_mode = 'strict_local_min';
% 
% scan = scan_sigma(problem, scan_opts);
% plot(scan.lamvec, scan.S)

%% body

% ---------------------------------
% 1. add project folders to path
% ---------------------------------

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(fullfile(project_root, 'cases'));
addpath(fullfile(project_root, 'polygon'));
addpath(fullfile(project_root, 'problem_builders'));
addpath(fullfile(project_root, 'core'));
addpath(fullfile(project_root, 'solvers'));
addpath(fullfile(project_root, 'minimizers'));
addpath(fullfile(project_root, 'plots'));

% ---------------------------------
% 2. build problem
% ---------------------------------
cfg = case_polygon_hshape();
cfg.qr_pivot = true;
% cfg = case_polygon_drum('left');
% cfg = case_polygon_lshape();
cfg.sampling.boundary_fun = @sample_boundary_chebyshev;
% TODO: spport @sample_interior_uniform (is it necessary??)
cfg.nb_per_edge = 50;
cfg.Mcorner = 50;
cfg.nI = 100;


cfg = case_polygon_drum('left');
cfg.qr_tau = 1e-13;      % no truncation, keep QB size fixed
cfg.nb_per_edge = 40;
cfg.Mcorner = 40;
cfg.nI = 50;
problem = build_polygon_problem(cfg);
sigma_fun = @(lam) problem.ops.sigma(lam);

% ---------------------------------
% A's numerical rank
% ---------------------------------
lam_list = 7:0.025:11;
A_op = problem.ops.A;
ss1 = [];
for lam=lam_list
    sig = svd(A_op(lam), "econ");

    ss1 = [ss1, sig(end)];
end
figure;
semilogy(lam_list, ss1(:), 'b-', 'LineWidth', 1.5);hold on;
% semilogy(lam_list, ss2(:)-ss1(:), 'b-', 'LineWidth', 1.5);
xlabel('\lambda');
ylabel('\sigma_{min}(A(\lambda))');
% title('smallest singular values of A(\lambda)');
% legend('smallest', 'Interpreter', 'tex');
grid on;
save_plot_eps('sigA_fullrank');
%%
% ---------------------------------
% 3. solver options
% ---------------------------------
opts = struct();

opts.scan = struct();
opts.scan.lamvec = 7:0.025:11;
opts.scan.detect_mode = 'strict_local_min';

opts.refine = struct();
opts.refine.bracket_halfwidth = 1;
opts.refine.sigma_cut = 5e-2;
opts.refine.minimizer = @minimizer_trisection;
opts.refine.minimizer_opts = struct( ...
    'tol_x', 1e-13, ...
    'tol_fun', 1e-12, ...
    'max_iter', 200, ...
    'max_fun_evals', 500, ...
    'n_pre', 10);

opts.refine.pre_refine_filter = [];
opts.refine.local_analyzer = [];

% ---------------------------------
% 4. solve
% ---------------------------------
result = solve_scan_refine(sigma_fun, opts);

% ---------------------------------
% 5. print summary
% ---------------------------------
print_candidate_summary(result);

% ---------------------------------
% 6. geometry plot
% ---------------------------------
geom = problem.data.geom;

figure;
plot(geom.V(:,1), geom.V(:,2), 'k-', 'LineWidth', 1.5); hold on;
plot(geom.B(:,1), geom.B(:,2), 'bo', 'MarkerSize', 4);
plot(geom.I(:,1), geom.I(:,2), 'r.', 'MarkerSize', 10);
axis equal;
grid on;
legend('polygon', 'boundary samples', 'interior samples', 'Location', 'best');
title(problem.name);
save_plot_eps('polygon_geometry');

% ---------------------------------
% 7. scan plot with refined candidates
% ---------------------------------
figure;
semilogy(result.scan.lamvec, result.scan.S, 'k-', 'LineWidth', 1.2); hold on;
grid on;
xlabel('\lambda');
% ylabel('\sigma(\lambda)');
% title('scan-refine result');

% % raw detected dips
% J = result.scan.candidate_idx;
% plot(result.scan.lamvec(J), result.scan.S(J), 'ko', 'MarkerFaceColor', 'y');

% accepted refined candidates
accepted = result.summary.accepted_mask;
cand_all = result.candidates;
cand_acc = cand_all(accepted);

if ~isempty(cand_acc)
    x = [cand_acc.refined_lambda];
    y = [cand_acc.refined_sigma];
    plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

    for i = 1:numel(cand_acc)
        xx =cand_acc(i).refined_lambda;
        yl = ylim;
        plot([xx xx], yl, 'r--');
    end
end

legend('scan', 'accepted refined', 'FontSize', 13, 'Location','best');
save_plot_eps('scan_refine_fullrank');
%%
lam = 12.335964909866;
QB = problem.ops.QB(12.336992698243);
[~,S,V] = svd(QB,"econ", "vector");
[~, idx] = min(S);
v = V(:, idx);
[n, m] = size(QB);
Q = QB*(eye(m) - v*v');
ss = svd(Q);
Q_op = @(lam) (problem.ops.QB(lam))*(eye(m) - v*v');
Q_op = @(lam) problem.ops.QB(lam);
lam_vec = 12.32:5e-5:12.34;
s1 = [];
s2 = [];
s3 = [];
for lam = lam_vec
    ss = svd(Q_op(lam),"econ");
    s1 = [s1, ss(end)];
    s2 = [s2, ss(end-1)];
    
    s3 = [s3, ss(end-2)];
end
figure;
semilogy(lam_vec, s1, 'k-', 'LineWidth', 1.5);hold on;
semilogy(lam_vec, s2, 'b-', 'LineWidth', 1.5);
semilogy(lam_vec, s3, 'g-', 'LineWidth', 1.5);
legend('\sigma_k','\sigma_{k-1}','\sigma_{k-2}','Location', 'best');
grid on;
save_plot_eps('singular_value_trace_cluster');
%%
lam0 = 12.337002090000;

h = 1e-4;

pair_opts = struct();
pair_opts.normalize_columns = cfg.normalize_columns;
pair_opts.qr_tau = cfg.qr_tau;
pair_opts.pivot = true;
pair_opts.sign_fix = true;

[dQB, QB_left, QB_right, pair_info] = approx_QB_derivative_from_Aop( ...
    problem.ops.A, problem.meta.mB, lam0, h, pair_opts);

fprintf('\npaired-Q example around lambda = %.12f\n', lam0);
fprintf('  rank at lambda-h          : %d\n', pair_info.rank_left);
fprintf('  rank at lambda+h          : %d\n', pair_info.rank_right);
fprintf('  common selected columns   : %d\n', pair_info.n_common);
fprintf('  gap at lambda-h           : %.3e\n', pair_info.gap_left);
fprintf('  gap at lambda+h           : %.3e\n', pair_info.gap_right);

if pair_info.n_common > 0
    fprintf('  ||QB(lambda+h)-QB(lambda-h)||_F / (2h) = %.3e\n', norm(dQB));
end
%%
a = opts.scan.lamvec(1);
b = opts.scan.lamvec(end);
accepted = result.summary.accepted_mask;
cand_all = result.candidates;
cand_acc = cand_all(accepted);

cluster_opts = struct();
cluster_opts.cluster_ratio = 100;
cluster_opts.fd_step = 1e-4;
cluster_opts.micro_max_points = 4001;

if isempty(cand_acc)
    fprintf('\nNo accepted refined candidates. Skip local cluster resolution.\n');
else
    fprintf('\n=== Local Cluster Resolution ===\n');
    for ic = 1:numel(cand_acc)
        cand = cand_acc(ic);
        report = resolve_local_cluster(problem, sigma_fun, cand, [a, b], opts.refine, cfg, cluster_opts);
        print_cluster_report(report);

        if report.cluster_detected && ~isempty(report.micro_scan.lamvec)
            figure;
            h_scan = plot(report.micro_scan.lamvec, report.micro_scan.S, 'k-', 'LineWidth', 1.2); hold on;
            grid on;
            h_line = xline(report.lam_star, 'b--', 'LineWidth', 1.0);
            h_base = plot(report.lam_star, report.sigma_min, 'bo', 'MarkerFaceColor', 'b');

            extra_mask = report.micro_summary.extra_mask;
            if any(extra_mask)
                extra_lam = [report.micro_summary.refined_candidates(extra_mask).refined_lambda];
                extra_sig = [report.micro_summary.refined_candidates(extra_mask).refined_sigma];
                h_extra = plot(extra_lam, extra_sig, 'ro', 'MarkerFaceColor', 'r');
            end
            
            set(gca, 'YScale', 'log')
            xlabel('\lambda');
            ylabel('\sigma_{min}(Q_B(\lambda))');
            title(sprintf('Local micro-scan near \\lambda_*'));
            if any(extra_mask)
                legend([h_scan, h_line, h_base, h_extra], ...
                    {'micro-scan', '\lambda_*', 'base dip', 'extra dips'}, ...
                    'Location', 'best');
            else
                legend([h_scan, h_line, h_base], ...
                    {'micro-scan', '\lambda_*', 'base dip'}, ...
                    'Location', 'best');
            end
            save_plot_eps(sprintf('local_micro_scan_%02d', ic));
        end
    end
end


%% ============================================================
% local function
% ============================================================
function s1 = smallest_sigma(Q_op, lam)
    s = svd(Q_op(lam), 'econ');
    s1 = s(end);
end


function report = resolve_local_cluster(problem, sigma_fun, cand, global_interval, refine_opts, cfg, cluster_opts)
    lam_star = cand.refined_lambda;
    L_macro = max(abs(cand.bracket - lam_star));

    [QB_star, factor_info] = build_reference_QB(problem, lam_star, cfg);
    [U, S] = svd(QB_star, 'econ');
    svals = diag(S);
    sigma_min = svals(end);

    k = find(svals < cluster_opts.cluster_ratio * sigma_min, 1, 'first');
    if isempty(k)
        k = numel(svals);
    end

    report = struct();
    report.lam_star = lam_star;
    report.sigma_min = sigma_min;
    report.cluster_index = k;
    report.cluster_detected = (k < numel(svals));
    report.n_small_singular = numel(svals) - k + 1;
    report.svals = svals;
    report.L_macro = L_macro;
    report.L_fine = [];
    report.fd_step = cluster_opts.fd_step;
    report.reference_rank = factor_info.rank;
    report.reference_cols = factor_info.selected_cols;
    report.micro_scan = struct('lamvec', [], 'S', [], 'candidate_idx', []);
    report.micro_summary = struct('refined_candidates', [], 'extra_mask', []);
    report.message = '';

    if ~report.cluster_detected
        report.message = 'No cluster detected at the refined lambda.';
        return;
    end

    QB_plus = build_QB_with_selected_cols(problem.ops.A, problem.meta.mB, lam_star + cluster_opts.fd_step, factor_info.selected_cols, cfg);
    delta_Q = QB_plus - QB_star;
    response = norm(delta_Q, 2);

    if response <= eps(class(response))
        report.message = 'Cluster detected, but the finite-difference response is numerically zero.';
        return;
    end

    L_fine = 0.5 * sigma_min * cluster_opts.fd_step / response;
    if ~isfinite(L_fine) || L_fine <= 0
        report.message = 'Cluster detected, but the adaptive local step is not finite.';
        return;
    end

    report.L_fine = L_fine;

    lam_left = max(global_interval(1), lam_star - L_macro);
    lam_right = min(global_interval(2), lam_star + L_macro);
    n_micro = max(3, ceil((lam_right - lam_left) / L_fine) + 1);

    if n_micro > cluster_opts.micro_max_points
        n_micro = cluster_opts.micro_max_points;
        report.message = sprintf('Micro-grid capped at %d points for efficiency.', n_micro);
    end

    lamvec = linspace(lam_left, lam_right, n_micro);
    S_micro = zeros(size(lamvec));
    for i = 1:numel(lamvec)
        S_micro(i) = sigma_fun(lamvec(i));
    end

    J = 2:numel(lamvec)-1;
    J = J(S_micro(J) < S_micro(J-1) & S_micro(J) < S_micro(J+1));

    micro_scan = struct();
    micro_scan.lamvec = lamvec;
    micro_scan.S = S_micro;
    micro_scan.candidate_idx = J;
    micro_scan.meta = struct();
    micro_scan.meta.detect_mode = 'strict_local_min';
    micro_scan.meta.n_scan_points = numel(lamvec);
    micro_scan.meta.n_candidates = numel(J);

    refined_candidates = refine_candidates(sigma_fun, micro_scan, refine_opts);
    extra_mask = local_extra_dip_mask(refined_candidates, lam_star, L_fine);

    report.micro_scan = micro_scan;
    report.micro_summary.refined_candidates = refined_candidates;
    report.micro_summary.extra_mask = extra_mask;

    if ~any(extra_mask)
        if isempty(report.message)
            report.message = 'Cluster detected, but no extra resolved dip was found in the micro-scan.';
        end
    elseif isempty(report.message)
        report.message = 'Cluster detected and at least one extra local dip was found.';
    end
end


function [QB, info] = build_reference_QB(problem, lam, cfg)
    A = problem.ops.A(lam);
    qr_opts = struct();
    qr_opts.normalize_columns = false;
    qr_opts.pivot = get_cfg_flag(cfg, 'qr_pivot');
    qr_opts.sign_fix = true;
    qr_opts.qr_tau = get_cfg_value(cfg, 'qr_tau', []);
    [QB, info] = build_QB_from_A(A, problem.meta.mB, qr_opts);
end


function QB = build_QB_with_selected_cols(A_op, mB, lam, selected_cols, cfg)
    A = A_op(lam);
    A = A(:, selected_cols);
    [Q, R] = qr(A, 0);

    d = sign(diag(R));
    d(d == 0) = 1;
    Q = Q * diag(d);

    QB = Q(1:mB, :);
end


function extra_mask = local_extra_dip_mask(candidates, lam_star, lam_tol)
    if isempty(candidates)
        extra_mask = false(0, 1);
        return;
    end

    extra_mask = false(numel(candidates), 1);
    for i = 1:numel(candidates)
        if ~candidates(i).accepted
            continue;
        end
        if abs(candidates(i).refined_lambda - lam_star) > lam_tol
            extra_mask(i) = true;
        end
    end
end


function print_cluster_report(report)
    fprintf('\n--- Cluster check at lambda* = %.15f ---\n', report.lam_star);
    fprintf('sigma_min              : %.6e\n', report.sigma_min);
    fprintf('cluster index k        : %d of %d\n', report.cluster_index, numel(report.svals));
    fprintf('small singular values  : %d\n', report.n_small_singular);
    fprintf('reference rank         : %d\n', report.reference_rank);

    if ~report.cluster_detected
        fprintf('status                 : no cluster\n');
        fprintf('message                : %s\n', report.message);
        return;
    end

    fprintf('status                 : cluster detected\n');
    fprintf('L_macro                : %.6e\n', report.L_macro);

    if isempty(report.L_fine)
        fprintf('message                : %s\n', report.message);
        return;
    end

    fprintf('L_fine                 : %.6e\n', report.L_fine);
    fprintf('micro points           : %d\n', numel(report.micro_scan.lamvec));

    n_extra = nnz(report.micro_summary.extra_mask);
    fprintf('extra dips found       : %d\n', n_extra);
    if n_extra == 0
        fprintf('extra dip status       : none resolved\n');
    else
        extra_cands = report.micro_summary.refined_candidates(report.micro_summary.extra_mask);
        for i = 1:numel(extra_cands)
            fprintf('  extra dip %d: lambda = %.15f, sigma = %.6e\n', ...
                i, extra_cands(i).refined_lambda, extra_cands(i).refined_sigma);
        end
    end

    if ~isempty(report.message)
        fprintf('message                : %s\n', report.message);
    end
end


function value = get_cfg_flag(cfg, name)
    value = isfield(cfg, name) && ~isempty(cfg.(name)) && cfg.(name);
end


function value = get_cfg_value(cfg, name, default_value)
    if isfield(cfg, name) && ~isempty(cfg.(name))
        value = cfg.(name);
    else
        value = default_value;
    end
end

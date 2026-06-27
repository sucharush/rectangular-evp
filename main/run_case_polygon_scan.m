clear; close all; clc;
rng(0);

%% body

% ---------------------------------
% 1. add project folders to path
% ---------------------------------

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(fullfile(project_root, 'cases'));
addpath(fullfile(project_root, 'polygon'));
addpath(fullfile(project_root, 'core'));
addpath(fullfile(project_root, 'solvers'));
addpath(fullfile(project_root, 'minimizers'));
addpath(fullfile(project_root, 'plots'));
% ---------------------------------
% 2. build problem
% ---------------------------------
cfg = case_polygon_hshape();
% cfg = case_polygon_lshape();
% cfg = case_polygon_drum('right');
cfg.qr_pivot = true;
cfg.qr_tau = 1e-13;      % scan truncation, matching the old scan logic
cfg.nb_per_edge = 60;
cfg.Mcorner = 60;
cfg.nI = 50;

do_plots = true;
do_singular_trace = false;
do_pair_derivative_example = true;

problem = build_polygon_problem(cfg);
sigma_fun_scan = @(lam) problem.ops.sigma(lam);

cfg_refine = cfg;
cfg_refine.qr_pivot = false;  % refinement: plain QR, no pivoting
cfg_refine.qr_tau = [];       % refinement: no rank truncation
sigma_fun_refine = @(lam) sigma_polygon(problem.data.geom, cfg_refine, lam);

% % ---------------------------------
% % A's numerical rank
% % ---------------------------------
% lam_list = 7:0.025:11;
% A_op = problem.ops.A;
% ss1 = [];
% for lam=lam_list
%     sig = svd(A_op(lam), "econ");
%
%     ss1 = [ss1, sig(end)];
% end
% figure;
% semilogy(lam_list, ss1(:), 'b-', 'LineWidth', 1.5);hold on;
% xlabel('\lambda');
% ylabel('\sigma_{min}(A(\lambda))');
% grid on;
% save_plot_eps('sigA_fullrank');
%%
% ---------------------------------
% 3. solver options
% ---------------------------------
opts = struct();

opts.scan = struct();
opts.scan.lamvec =  1:0.025:6;
opts.scan.lamvec =  10:0.025:13;
opts.scan.detect_mode = 'strict_local_min';

opts.refine = struct();
opts.refine.bracket_halfwidth = 1;
opts.refine.sigma_cut = 1e-3;
opts.refine.minimizer = @minimizer_golden_section;
% opts.refine.minimizer = @minimizer_fminsearch;
opts.refine.minimizer_opts = struct( ...
    'tol_x', 1e-13, ...
    'tol_fun', 1e-12, ...
    'max_iter', 200, ...
    'max_fun_evals', 500, ...
    'n_pre', 0);

% opts.refine.minimizer = @minimizer_aaa_real;
% 
% opts.refine.minimizer = @minimizer_aaa_ellipse;
% 
% opts.refine.minimizer_opts = struct( ...
%     'nZ', 100, ...
%     'delta', 1e-13, ...
%     'imag_tol', 1e-3, ...
%     'mmax', 100, ...
%     'eval_tol', 1e-13, ...
%     'rho', 1.05);

opts.refine.pre_refine_filter = [];
opts.refine.local_analyzer = [];

% ---------------------------------
% 4. solve
% ---------------------------------
scan = scan_sigma(sigma_fun_scan, opts.scan);
candidates = refine_candidates(sigma_fun_refine, scan, opts.refine);
summary = summarize_candidates(candidates);

result = struct();
result.method = 'scan_refine_old_logic';
result.scan = scan;
result.candidates = candidates;
result.summary = summary;

% % save("result_left.mat");
% load("result_left.mat");
% result = result_left;

% ---------------------------------
% 5. print summary
% ---------------------------------
print_candidate_summary(result);

% ---------------------------------
% 6. plots
% ---------------------------------
if do_plots
    % plot_polygon_geometry(problem.data.geom, ...
    %     'plot_title', problem.name, ...
    %     'save_name', 'polygon_geometry');
    % scan-refine plot is drawn at the end so it reflects any cluster merge
end
result_left = result;
% %%
% % ---------------------------------
% % 2. build problem
% % ---------------------------------
% cfg = case_polygon_hshape();
% cfg = case_polygon_lshape();
% cfg = case_polygon_drum('right');
% cfg.qr_pivot = true;
% cfg.qr_tau = 1e-13;      % scan truncation, matching the old scan logic
% cfg.nb_per_edge = 140;
% cfg.Mcorner = 140;
% cfg.nI = 50;
%
% do_plots = true;
% do_singular_trace = false;
% do_pair_derivative_example = false;
%
% problem = build_polygon_problem(cfg);
% sigma_fun_scan = @(lam) problem.ops.sigma(lam);
%
% cfg_refine = cfg;
% cfg_refine.qr_tau = [];  % old refine logic: QR without truncation
% sigma_fun_refine = @(lam) sigma_polygon(problem.data.geom, cfg_refine, lam);
%
% % % ---------------------------------
% % % A's numerical rank
% % % ---------------------------------
% % lam_list = 7:0.025:11;
% % A_op = problem.ops.A;
% % ss1 = [];
% % for lam=lam_list
% %     sig = svd(A_op(lam), "econ");
% %
% %     ss1 = [ss1, sig(end)];
% % end
% % figure;
% % semilogy(lam_list, ss1(:), 'b-', 'LineWidth', 1.5);hold on;
% % xlabel('\lambda');
% % ylabel('\sigma_{min}(A(\lambda))');
% % grid on;
% % save_plot_eps('sigA_fullrank');
% %%
% % ---------------------------------
% % 3. solver options
% % ---------------------------------
% opts = struct();
%
% opts.scan = struct();
% opts.scan.lamvec = 1:0.025:6;
% opts.scan.detect_mode = 'strict_local_min';
%
% opts.refine = struct();
% opts.refine.bracket_halfwidth = 1;
% opts.refine.sigma_cut = 1e-3;
% opts.refine.minimizer = @minimizer_fminsearch;
% opts.refine.minimizer = @minimizer_golden_section;
% opts.refine.minimizer_opts = struct( ...
%     'tol_x', 1e-13, ...
%     'tol_fun', 1e-12, ...
%     'max_iter', 200, ...
%     'max_fun_evals', 500, ...
%     'n_pre', 0);
%
%
% opts.refine.pre_refine_filter = [];
% opts.refine.local_analyzer = [];
%
% % ---------------------------------
% % 4. solve
% % ---------------------------------
% scan = scan_sigma(sigma_fun_scan, opts.scan);
% candidates = refine_candidates(sigma_fun_refine, scan, opts.refine);
% summary = summarize_candidates(candidates);
%
% result = struct();
% result.method = 'scan_refine_old_logic';
% result.scan = scan;
% result.candidates = candidates;
% result.summary = summary;
%
% % % save("result_left.mat");
% % load("result_left.mat");
% % result = result_left;
%
% % ---------------------------------
% % 5. print summary
% % ---------------------------------
% print_candidate_summary(result);
% result_right = result;
%
% %%
%
% % To compare left/right scans:
% % 1. Run this script for case_polygon_drum('left'), then set result_left = result.
% % 2. Run it again for case_polygon_drum('right'), then set result_right = result.
% % 3. Uncomment:
% figure;
% h_left = semilogy(result_left.scan.lamvec, result_left.scan.S, 'k-', 'LineWidth', 1.2); hold on;
% h_right = semilogy(result_right.scan.lamvec, result_right.scan.S, 'b-', 'LineWidth', 1.2);
% h_leg = [h_left, h_right];
% leg_names = {'GWW1 scan', 'GWW2 scan'};
%
% if ~isempty(result_left.summary.eigs)
%     h_left_ref = semilogy(result_left.summary.eigs, result_left.summary.sigmins, ...
%         'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
%     h_leg(end+1) = h_left_ref;
%     leg_names{end+1} = 'GWW1 refined';
% end
%
% if ~isempty(result_right.summary.eigs)
%     h_right_ref = semilogy(result_right.summary.eigs, result_right.summary.sigmins, ...
%         'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
%     h_leg(end+1) = h_right_ref;
%     leg_names{end+1} = 'GWW2 refined';
% end
%
% yl = ylim;
% for i = 1:numel(result_left.summary.eigs)
%     plot([result_left.summary.eigs(i), result_left.summary.eigs(i)], yl, 'k--', 'LineWidth', 1.0);
% end
% for i = 1:numel(result_right.summary.eigs)
%     plot([result_right.summary.eigs(i), result_right.summary.eigs(i)], yl, 'b--', 'LineWidth', 1.0);
% end
% ylim(yl);
%
% grid on;
% xlabel('\lambda');
% ylabel('\sigma_{min}(Q_B(\lambda))');
% legend(h_leg, leg_names, 'FontSize', 13, 'Location', 'best');
% save_plot_eps('polygon_drum_left_right_scan');
%%
if do_singular_trace
    Q_op = @(lam) problem.ops.QB(lam);
    lam_vec = 12.32:5e-5:12.34;
    n_trace = 3;
    sigma_trace = zeros(numel(lam_vec), n_trace);
    for ilam = 1:numel(lam_vec)
        svals = svd(Q_op(lam_vec(ilam)), "econ");
        sigma_trace(ilam, :) = svals(end:-1:end-n_trace+1).';
    end
    plot_singular_value_trace(lam_vec, sigma_trace, ...
        'labels', {'\sigma_k', '\sigma_{k-1}', '\sigma_{k-2}'}, ...
        'save_name', 'singular_value_trace_cluster');
end
%%
if do_pair_derivative_example
    lam0 = 12.337050160494;
    h = 1e-4;

    pair_opts = struct();
    pair_opts.normalize_columns = cfg.normalize_columns;
    pair_opts.qr_tau = cfg.qr_tau;
    pair_opts.pivot = true;
    pair_opts.sign_fix = true;

    [dQB, ~, ~, pair_info] = approx_QB_derivative_from_Aop( ...
        problem.ops.A, problem.meta.mB, lam0, h, pair_opts);

    fprintf('\npaired-Q example around lambda = %.12f\n', lam0);
    fprintf('  rank at lambda-h          : %d\n', pair_info.rank_left);
    fprintf('  rank at lambda+h          : %d\n', pair_info.rank_right);
    fprintf('  common selected columns   : %d\n', pair_info.n_common);
    fprintf('  gap at lambda-h           : %.3e\n', pair_info.gap_left);
    fprintf('  gap at lambda+h           : %.3e\n', pair_info.gap_right);

    if pair_info.n_common > 0
        fprintf('  ||QB(lambda+h)-QB(lambda-h)||_2 / (2h) = %.3e\n', norm(dQB, 2));
    end
end
%% post-processing
a = opts.scan.lamvec(1);
b = opts.scan.lamvec(end);
accepted = result.summary.accepted_mask;
acc_idx = find(accepted);
cand_acc = result.candidates(accepted);

cluster_opts = struct();
cluster_opts.cluster_ratio = 100;
cluster_opts.fd_step = 1e-4;
cluster_opts.micro_max_points = 4001;
cluster_opts.dedup_factor = 8;   % dedup tol = dedup_factor * L_fine

if isempty(cand_acc)
    fprintf('\nNo accepted refined candidates. Skip local cluster resolution.\n');
else
    fprintf('\n=== Local Cluster Resolution ===\n');
    any_replaced = false;
    for ic = 1:numel(cand_acc)
        cand = cand_acc(ic);
        report = resolve_local_cluster(problem, sigma_fun_scan, sigma_fun_refine, ...
            cand, [a, b], opts.refine, cfg, cluster_opts);
        print_cluster_report(report);

        % If a deeper re-detection merged with this dip, replace the candidate.
        if report.merged_replaced
            gi = acc_idx(ic);
            result.candidates(gi).refined_lambda = report.merged_lambda;
            result.candidates(gi).refined_sigma = report.merged_sigma;
            any_replaced = true;
        end

        if do_plots && report.cluster_detected && ~isempty(report.micro_scan.lamvec)
            plot_local_micro_scan_report(report, ...
                'save_name', sprintf('local_micro_scan_%02d', ic));
        end
    end

    if any_replaced
        result.summary = summarize_candidates(result.candidates);
        fprintf('\n=== Candidate summary after cluster merge: ===\n');
        print_candidate_summary(result);
    end
end

% scan-refine plot, drawn after post-processing so it reflects cluster merges
if do_plots
    plot_scan_refine_result(result, ...
        'problem_name', problem.name, ...
        'save_name', 'raw_cluster');
end


%% ============================================================
% local function
% ============================================================
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
    fprintf('response ||dQB||_2     : %.6e\n', report.response);

    if isempty(report.L_fine)
        fprintf('message                : %s\n', report.message);
        return;
    end

    fprintf('L_fine                 : %.6e\n', report.L_fine);
    fprintf('dedup tol              : %.6e\n', report.dedup_tol);
    fprintf('micro points           : %d\n', numel(report.micro_scan.lamvec));

    if report.merged_replaced
        fprintf('merged dip             : lambda = %.15f, sigma = %.6e (replaced original sigma %.6e)\n', ...
            report.merged_lambda, report.merged_sigma, report.refined_sigma_original);
    else
        fprintf('merged dip             : lambda = %.15f, sigma = %.6e (original retained, deepest of %d merged)\n', ...
            report.merged_lambda, report.merged_sigma, report.n_merged);
    end

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

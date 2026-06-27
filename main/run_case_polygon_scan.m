clear; close all; clc;
rng(0);

% ============================================================
% Pick the experiment (see experiments/polygon_experiment_specs.m):
%   lshape_scan
%   gww1_scan
%   gww2_scan
%   hshape_width_0p3_scan
%   hshape_width_0p08_cluster   (scan + cluster-aware post-processing)
%   gww_pair_scan               (GWW1 vs GWW2 overlay)
% ============================================================
spec_name = 'hshape_width_0p08_cluster';

% ---------------------------------
% project paths
% ---------------------------------
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(this_dir);

addpath(fullfile(project_root, 'cases'));
addpath(fullfile(project_root, 'polygon'));
addpath(fullfile(project_root, 'core'));
addpath(fullfile(project_root, 'solvers'));
addpath(fullfile(project_root, 'minimizers'));
addpath(fullfile(project_root, 'plots'));
addpath(fullfile(project_root, 'experiments'));

do_plots = true;

% ---------------------------------
% shared options (same for every experiment)
% QR policy invariant: scan = CPQR + qr_tau = 1e-13; refine = plain QR + qr_tau = []
% ---------------------------------
refine_opts = struct();
refine_opts.bracket_halfwidth = 1;
refine_opts.sigma_cut = 1e-3;
refine_opts.minimizer = @minimizer_golden_section;
refine_opts.minimizer_opts = struct( ...
    'tol_x', 1e-13, ...
    'tol_fun', 1e-12, ...
    'max_iter', 200, ...
    'max_fun_evals', 500, ...
    'n_pre', 0);
refine_opts.pre_refine_filter = [];
refine_opts.local_analyzer = [];

cluster_opts = struct();
cluster_opts.cluster_ratio = 100;
cluster_opts.fd_step = 1e-4;
cluster_opts.micro_max_points = 4001;
cluster_opts.dedup_factor = 8;   % dedup tol = dedup_factor * L_fine

% ---------------------------------
% run the selected experiment
% ---------------------------------
spec = polygon_experiment_specs(spec_name);

if isfield(spec, 'pair_members')
    results = cell(numel(spec.pair_members), 1);
    for i = 1:numel(spec.pair_members)
        member = polygon_experiment_specs(spec.pair_members{i});
        results{i} = run_scan_experiment_one(member, refine_opts, cluster_opts, false);
    end
    if do_plots
        plot_pair_overlay(results, spec);
    end
else
    result = run_scan_experiment_one(spec, refine_opts, cluster_opts, do_plots);
end


%% ============================================================
% local functions
% ============================================================
function result = run_scan_experiment_one(spec, refine_opts, cluster_opts, do_plots)
    % geometry cfg + scan QR policy
    cfg = spec.case_fun();
    cfg.qr_pivot = true;
    cfg.qr_tau = 1e-13;
    cfg.Mcorner = spec.Mcorner;
    cfg.nb_per_edge = spec.nb_per_edge;
    cfg.nI = spec.nI;

    problem = build_polygon_problem(cfg);
    sigma_fun_scan = @(lam) problem.ops.sigma(lam);

    % refine QR policy: plain QR, no truncation
    cfg_refine = cfg;
    cfg_refine.qr_pivot = false;
    cfg_refine.qr_tau = [];
    sigma_fun_refine = @(lam) sigma_polygon(problem.data.geom, cfg_refine, lam);

    scan_opts = struct();
    scan_opts.lamvec = spec.interval(1):spec.step:spec.interval(2);
    scan_opts.detect_mode = 'strict_local_min';

    scan = scan_sigma(sigma_fun_scan, scan_opts);
    candidates = refine_candidates(sigma_fun_refine, scan, refine_opts);
    summary = summarize_candidates(candidates);

    result = struct();
    result.spec_name = spec.name;
    result.scan = scan;
    result.candidates = candidates;
    result.summary = summary;

    fprintf('\n### experiment: %s ###\n', spec.name);
    print_candidate_summary(result);

    % cluster-aware post-processing (only the thin H-shape)
    if spec.cluster_enable
        result = resolve_clusters(result, problem, sigma_fun_scan, sigma_fun_refine, ...
            cfg, refine_opts, cluster_opts, do_plots);
    end

    if do_plots
        plot_scan_refine_result(result, ...
            'problem_name', problem.name, ...
            'save_name', spec.plot_name);
    end
end


function result = resolve_clusters(result, problem, sigma_fun_scan, sigma_fun_refine, ...
        cfg, refine_opts, cluster_opts, do_plots)
    a = result.scan.lamvec(1);
    b = result.scan.lamvec(end);
    accepted = result.summary.accepted_mask;
    acc_idx = find(accepted);
    cand_acc = result.candidates(accepted);

    if isempty(cand_acc)
        fprintf('\nNo accepted refined candidates. Skip local cluster resolution.\n');
        return;
    end

    fprintf('\n=== Local Cluster Resolution ===\n');
    any_replaced = false;
    for ic = 1:numel(cand_acc)
        cand = cand_acc(ic);
        report = resolve_local_cluster(problem, sigma_fun_scan, sigma_fun_refine, ...
            cand, [a, b], refine_opts, cfg, cluster_opts);
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


function plot_pair_overlay(results, spec)
    colors = {'k', 'b'};
    figure;
    h_leg = gobjects(0);
    leg_names = {};

    for i = 1:numel(results)
        res = results{i};
        c = colors{mod(i-1, numel(colors)) + 1};

        h = semilogy(res.scan.lamvec, res.scan.S, [c '-'], 'LineWidth', 1.2);
        hold on;
        h_leg(end+1) = h; %#ok<AGROW>
        leg_names{end+1} = sprintf('%s scan', spec.pair_labels{i}); %#ok<AGROW>

        if ~isempty(res.summary.eigs)
            hr = semilogy(res.summary.eigs, res.summary.sigmins, [c 'o'], ...
                'MarkerFaceColor', c, 'MarkerSize', 5);
            h_leg(end+1) = hr; %#ok<AGROW>
            leg_names{end+1} = sprintf('%s refined', spec.pair_labels{i}); %#ok<AGROW>
        end
    end

    yl = ylim;
    for i = 1:numel(results)
        c = colors{mod(i-1, numel(colors)) + 1};
        for e = results{i}.summary.eigs(:)'
            plot([e e], yl, [c '--'], 'LineWidth', 1.0);
        end
    end
    ylim(yl);

    grid on;
    xlabel('\lambda');
    ylabel('\sigma_{min}(Q_B(\lambda))');
    legend(h_leg, leg_names, 'FontSize', 13, 'Location', 'best');
    save_plot_eps(spec.plot_name);
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

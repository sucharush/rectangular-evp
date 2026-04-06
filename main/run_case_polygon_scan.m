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
% cfg = case_polygon_hshape();
cfg = case_polygon_lshape();
cfg.sampling.boundary_fun = @sample_boundary_chebyshev;
% TODO: spport @sample_interior_uniform (is it necessary??)
cfg.nb_per_edge = 140;
cfg.Mcorner = 140;
cfg.nI = 50;
problem = build_polygon_problem(cfg);
sigma_fun = @(lam) problem.ops.sigma(lam);

% ---------------------------------
% 3. solver options
% ---------------------------------
opts = struct();

opts.scan = struct();
opts.scan.lamvec = 4:0.01:13;
opts.scan.detect_mode = 'strict_local_min';

opts.refine = struct();
opts.refine.bracket_halfwidth = 1;
opts.refine.sigma_cut = 1e-2;
opts.refine.minimizer = @minimizer_hybrid;
opts.refine.minimizer_opts = struct( ...
    'tol_x', 1e-13, ...
    'tol_fun', 1e-12, ...
    'max_iter', 200, ...
    'max_fun_evals', 500, ...
    'n_pre', 5);

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

% ---------------------------------
% 7. scan plot with refined candidates
% ---------------------------------
figure;
plot(result.scan.lamvec, result.scan.S, 'k-', 'LineWidth', 1.2); hold on;
grid on;
xlabel('\lambda');
ylabel('\sigma(\lambda)');
title('scan-refine result');

% raw detected dips
J = result.scan.candidate_idx;
plot(result.scan.lamvec(J), result.scan.S(J), 'ko', 'MarkerFaceColor', 'y');

% accepted refined candidates
accepted = result.summary.accepted_mask;
cand_all = result.candidates;
cand_acc = cand_all(accepted);

if ~isempty(cand_acc)
    x = [cand_acc.refined_lambda];
    y = [cand_acc.refined_sigma];
    plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

    for i = 1:numel(cand_acc)
        xx = cand_acc(i).refined_lambda;
        yl = ylim;
        plot([xx xx], yl, 'r--');
    end
end

legend('scan', 'raw dips', 'accepted refined', 'Location', 'best');

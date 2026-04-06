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

cfg = case_polygon_lshape();
cfg.qr_tau = 0; % no column pivoting to avoid different QB size
cfg.nb_per_edge = 15;
cfg.Mcorner = 15;
cfg.nI = 50;
%----------------------------------------
% TODO: why does col pivoting change QB's dim???
%----------------------------------------
problem = build_polygon_problem(cfg);

%----------------------------------------
% try aaa on A(lam): don't forget to change it back....
problem.ops.F = problem.ops.A;
%----------------------------------------

F_op = problem.ops.F;

opts = struct();
opts.interval = [2, 4];

opts.aaa = struct();
opts.aaa.nCand = 500;
opts.aaa.aaa_tol = 1e-12;
opts.aaa.mmax = 80;
opts.aaa.seed = 1;
opts.aaa.ell = 4;
opts.aaa.max_norm = true;
opts.aaa.svd_update = true;
opts.aaa.beta_mode = 'ones';
opts.aaa.rho = 1.05;
opts.aaa.sample_mode = 'ellipse';

opts.verify = struct();
opts.verify.do_scan = true;
opts.verify.lam_grid = linspace(2, 4, 100);

result = solve_aaa(F_op, opts);

disp(result.logs)

figure;
plot(result.verification.lam_grid, result.verification.sigma_original, 'k-', 'LineWidth', 1.2); hold on;
plot(result.verification.lam_grid, result.verification.sigma_bary, '--', 'LineWidth', 1.2);
plot(result.verification.lam_grid, result.verification.sigma_newton, ':', 'LineWidth', 1.2);
grid on;
xlabel('\lambda');
ylabel('tracked singular value');
% title(sprintf('AAA verification: %s', result.problem_name), 'Interpreter', 'none');
legend('original', 'barycentric', 'newton', 'Location', 'best');
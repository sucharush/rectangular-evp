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
% cfg = case_polygon_lshape();
cfg.sampling.boundary_fun = @sample_boundary_chebyshev;
% TODO: spport @sample_interior_uniform (is it necessary??)
cfg.nb_per_edge = 40;
cfg.Mcorner = 40;
cfg.nI = 50;
problem = build_polygon_problem(cfg);
sigma_fun = @(lam) problem.ops.sigma(lam);

% ---------------------------------
% 3. solver options
% ---------------------------------
opts = struct();

opts.scan = struct();
opts.scan.lamvec = 4:0.05:13;
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
%%
lam = 12.335964909866;
QB = problem.ops.QB(12.335964909866);
[~,S,V] = svd(QB,"econ", "vector");
[~, idx] = min(S);
v = V(:, idx);
[n, m] = size(QB);
Q = QB*(eye(m) - v*v');
ss = svd(Q);
Q_op = @(lam) (problem.ops.QB(lam))*(eye(m) - v*v');
Q_op = @(lam) problem.ops.QB(lam);
lam_vec = 12.3:1e-4:12.4;
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
%%
a = 12.3;
b = 12.4;
lam0 = 12.336987585210;

Q_op = @(lam) problem.ops.QB(lam);

abs_tol = 1e-2;
close_ratio = 10;
n_init = 21;
max_refine = 6;
min_width = 1e-8;

kmax_keep = 8;              % how many smallest singular values to inspect
local_half_width = 0.01;    % local window around lam0 for tracing
sep_tol_lam = 1e-5;         % tolerance for saying two traced minima coincide
dip_tol = abs_tol;          % threshold for accepting an extra sigma_min dip

sigma_min_fun = @(lam) smallest_sigma(Q_op, lam);

%% ============================================================
% Step 0. singular values at lam0
% ============================================================
s0 = svd(Q_op(lam0), 'econ');
kmax = min(kmax_keep, numel(s0));
s_small0 = flipud(s0(end-kmax+1:end));   % ascending

disp('small singular values at lam0 = ');
disp(s_small0(:).');

idx_abs = find(s_small0 <= abs_tol);
if isempty(idx_abs)
    error('No singular values below abs_tol at lam0.');
end

fprintf('absolute-small indices at lam0 = ');
disp(idx_abs);

%% ============================================================
% Step 1. grouping at lam0
% ============================================================
groups0 = {};
g = 1;
groups0{g} = idx_abs(1);

for t = 2:numel(idx_abs)
    j_prev = idx_abs(t-1);
    j_curr = idx_abs(t);
    ratio_tc = s_small0(j_curr) / max(s_small0(j_prev), eps);

    if ratio_tc <= close_ratio
        groups0{g}(end+1) = j_curr;
    else
        g = g + 1;
        groups0{g} = j_curr;
    end
end

fprintf('groups at lam0:\n');
for g = 1:numel(groups0)
    fprintf('  Group %d: ', g);
    fprintf('%d ', groups0{g});
    fprintf('\n');
end

active_group = groups0{1};
other_groups = groups0(2:end);
r_keep = idx_abs(end);

fprintf('active group = ');
disp(active_group);

%% ============================================================
% Step 2. trace only the active group locally near lam0
%         and judge separability of the known dip
% ============================================================
a_loc = max(a, lam0 - local_half_width);
b_loc = min(b, lam0 + local_half_width);

lam_active = unique(sort([linspace(a_loc, b_loc, n_init), lam0]));

for depth = 1:max_refine
    n = numel(lam_active);
    S = NaN(n, r_keep);

    for i = 1:n
        s = svd(Q_op(lam_active(i)), 'econ');
        s = flipud(s(end-r_keep+1:end));
        S(i,:) = s(:).';
    end

    new_pts = [];

    for i = 1:n-1
        l = lam_active(i);
        r = lam_active(i+1);
        if (r-l) < min_width
            continue;
        end

        m = 0.5*(l+r);
        if any(abs(lam_active - m) < 1e-15)
            continue;
        end

        s_mid = svd(Q_op(m), 'econ');
        s_mid = flipud(s_mid(end-r_keep+1:end));

        refine = false;

        % refine if any active-group branch shows a valley
        for j = active_group
            if (s_mid(j) < S(i,j)) && (s_mid(j) < S(i+1,j))
                refine = true;
                break;
            end
        end

        % also refine if active-group ordering looks unstable
        if ~refine && numel(active_group) >= 2
            left_order = S(i,active_group);
            right_order = S(i+1,active_group);
            mid_order = s_mid(active_group);

            if any(diff(left_order).*diff(right_order) < 0) || ...
               any(diff(left_order).*diff(mid_order) < 0)
                refine = true;
            end
        end

        if refine
            new_pts(end+1) = m; %#ok<AGROW>
        end
    end

    if isempty(new_pts)
        break;
    end

    lam_active = unique(sort([lam_active, new_pts]));
end

n = numel(lam_active);
S_active = NaN(n, r_keep);
for i = 1:n
    s = svd(Q_op(lam_active(i)), 'econ');
    s = flipud(s(end-r_keep+1:end));
    S_active(i,:) = s(:).';
end

figure;
hold on;
for j = active_group
    plot(lam_active, S_active(:,j), 'o-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('\\sigma_%d', j));
end
xline(lam0, '--k', 'HandleVisibility', 'off');
xlabel('\lambda');
ylabel('singular value');
title('Active group near known dip');
legend('Location', 'best');
grid on;

fprintf('\n=== Active group diagnostics ===\n');
lam_min_active = zeros(numel(active_group),1);
sig_min_active = zeros(numel(active_group),1);

for p = 1:numel(active_group)
    j = active_group(p);
    [sig_min_active(p), idxm] = min(S_active(:,j));
    lam_min_active(p) = lam_active(idxm);
    fprintf('sigma_%d min at %.15f, value %.6e\n', ...
        j, lam_min_active(p), sig_min_active(p));
end

if numel(active_group) == 1
    fprintf('Interpretation: active group is simple.\n');
else
    lam_spread = max(lam_min_active) - min(lam_min_active);
    sig_spread = max(sig_min_active) - min(sig_min_active);

    fprintf('spread of traced minimizers in lambda = %.6e\n', lam_spread);
    fprintf('spread of minimum values inside active group = %.6e\n', sig_spread);

    if lam_spread <= sep_tol_lam
        fprintf('Interpretation: active group behaves like one unresolved local cluster.\n');
    else
        fprintf('Interpretation: active group appears separable.\n');
    end
end

%% ============================================================
% Step 3. use other groups only to mark suspicious intervals
%         where sigma_min may have another dip
% ============================================================
lam_scan = linspace(a, b, n_init);
S_scan = NaN(numel(lam_scan), r_keep);

for i = 1:numel(lam_scan)
    s = svd(Q_op(lam_scan(i)), 'econ');
    s = flipud(s(end-r_keep+1:end));
    S_scan(i,:) = s(:).';
end

for depth = 1:max_refine
    new_pts = [];

    for i = 1:numel(lam_scan)-1
        l = lam_scan(i);
        r = lam_scan(i+1);
        if (r-l) < min_width
            continue;
        end

        m = 0.5*(l+r);
        if any(abs(lam_scan - m) < 1e-15)
            continue;
        end

        s_mid = svd(Q_op(m), 'econ');
        s_mid = flipud(s_mid(end-r_keep+1:end));

        refine = false;

        % ordinary sigma_min valley detection
        if (s_mid(1) < S_scan(i,1)) && (s_mid(1) < S_scan(i+1,1))
            refine = true;
        end

        % use other groups only as detectors of suspicious regions
        if ~refine
            for gg = 1:numel(other_groups)
                G = other_groups{gg};
                for j = G
                    if (s_mid(j) < S_scan(i,j)) && (s_mid(j) < S_scan(i+1,j))
                        refine = true;
                        break;
                    end
                end
                if refine
                    break;
                end
            end
        end

        if refine
            new_pts(end+1) = m; %#ok<AGROW>
        end
    end

    if isempty(new_pts)
        break;
    end

    lam_scan = unique(sort([lam_scan, new_pts]));
    S_scan = NaN(numel(lam_scan), r_keep);

    for i = 1:numel(lam_scan)
        s = svd(Q_op(lam_scan(i)), 'econ');
        s = flipud(s(end-r_keep+1:end));
        S_scan(i,:) = s(:).';
    end
end

sig1 = S_scan(:,1);

figure;
hold on;
plot(lam_scan, sig1, 'o-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{min}');
xline(lam0, '--k', 'known dip', 'HandleVisibility', 'off');
xlabel('\lambda');
ylabel('\sigma_{min}(Q_B(\lambda))');
title('Refined scan for possible additional dips of \sigma_{min}');
legend('Location', 'best');
grid on;

%% ============================================================
% Step 4. define suspicious intervals and solve sigma_min there
%         IMPORTANT: final lambda must minimize sigma_min
% ============================================================
idx_susp = find(sig1(2:end-1) <= sig1(1:end-2) & sig1(2:end-1) <= sig1(3:end)) + 1;

cand_intervals = [];
for t = 1:numel(idx_susp)
    k = idx_susp(t);
    l = lam_scan(max(k-1,1));
    r = lam_scan(min(k+1,numel(lam_scan)));

    % ignore the known dip near lam0
    if (l <= lam0) && (lam0 <= r)
        continue;
    end

    cand_intervals(end+1,:) = [l, r]; %#ok<AGROW>
end

fprintf('\n=== suspicious intervals for extra sigma_min dips ===\n');
disp(cand_intervals);

extra_dips = [];

for t = 1:size(cand_intervals,1)
    l = cand_intervals(t,1);
    r = cand_intervals(t,2);

    [lam_star, sig_star] = fminbnd(sigma_min_fun, l, r);

    fprintf('interval [% .15f, % .15f] -> lam_star = %.15f, sigma_min = %.6e\n', ...
        l, r, lam_star, sig_star);

    if sig_star <= dip_tol
        extra_dips(end+1,:) = [lam_star, sig_star, l, r]; %#ok<AGROW>
    end
end

fprintf('\n=== accepted extra dips from sigma_min minimization ===\n');
disp(extra_dips);

%% ============================================================
% Step 5. inspect each accepted extra dip locally
%         if small sigmas cluster there, report it
% ============================================================
for kk = 1:size(extra_dips,1)
    lamk = extra_dips(kk,1);

    sk = svd(Q_op(lamk), 'econ');
    kmax_k = min(kmax_keep, numel(sk));
    s_smallk = flipud(sk(end-kmax_k+1:end));

    idx_abs_k = find(s_smallk <= abs_tol);

    fprintf('\n=== local inspection at extra dip %d ===\n', kk);
    fprintf('lam = %.15f\n', lamk);
    fprintf('small singular values = ');
    disp(s_smallk(:).');

    if isempty(idx_abs_k)
        fprintf('No absolute-small singular values found under abs_tol.\n');
        continue;
    end

    groups_k = {};
    g = 1;
    groups_k{g} = idx_abs_k(1);

    for t = 2:numel(idx_abs_k)
        jp = idx_abs_k(t-1);
        jc = idx_abs_k(t);
        ratio_tc = s_smallk(jc) / max(s_smallk(jp), eps);

        if ratio_tc <= close_ratio
            groups_k{g}(end+1) = jc;
        else
            g = g + 1;
            groups_k{g} = jc;
        end
    end

    fprintf('groups at this extra dip:\n');
    for g = 1:numel(groups_k)
        fprintf('  Group %d: ', g);
        fprintf('%d ', groups_k{g});
        fprintf('\n');
    end

    % optional: if first group has size > 1, do another local trace there
    if numel(groups_k{1}) > 1
        fprintf('This extra dip also has a clustered active group; local tracing is recommended.\n');
    end
end

%% ============================================================
% local function
% ============================================================
function s1 = smallest_sigma(Q_op, lam)
    s = svd(Q_op(lam), 'econ');
    s1 = s(end);
end
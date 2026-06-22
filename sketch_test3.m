% clear; close all; clc;
rng(0);

%% Setup
n0 = 228;
vals = [-n0:-0.1, 0.1:0.1:2.8];
n = numel(vals);

[Q0, ~] = qr(randn(n));
A = Q0 * diag(vals) * Q0';

B = 5e-3 * randn(n);
B = B + B';
B = Q0 * B * Q0';

AA = A;
lam_true = 0.1;

iiMax_list = 1:20;
plot_dims = [10, 20];
n_rep = 30;
sketch_type = 'gaussian'; % 'gaussian', 'haar', 'srht'
% sketch_type = 'srht';
sketch_type = lower(char(sketch_type));
if ~ismember(sketch_type, {'gaussian', 'haar', 'srht'})
    error('Unknown sketch_type "%s". Use gaussian, haar, or srht.', sketch_type);
end

rad = 0.2;
h = 1e-3;
lam_grid = (0.2-rad:h:rad).';


[VAA, evalsAA] = eig(AA, 'vector');
[~, ind_true] = min(abs(evalsAA));
v_true = VAA(:, ind_true);
v_true = v_true / norm(v_true);

iiMax_target = iiMax_list(end);
Ufull = zeros(n, iiMax_target);
for ii = 1:iiMax_target
    [V, evals] = eig(A + ii * B, 'vector');
    [~, ind] = min(abs(evals));
    Ufull(:, ii) = V(:, ind);
end
subspace(Ufull, v_true)


ell_labels = {'l=r', 'l=2r', 'l=4r', 'l=n/2', 'l=n-r', 'l=n-1', 'l=n'};
% ell_labels = {'l=r', 'l=n/2', 'l=n-4r','l=n-2r','l=n-r', 'l=n-1', 'l=n'};
n_cases = numel(iiMax_list);
n_specs = numel(ell_labels);

dimQ = nan(n_cases, 1);
sub_err = nan(n_cases, 1);
abserr_full = nan(n_cases, 1);
angle_full = nan(n_cases, 1);
spec_valid = false(n_cases, n_specs);
abserr_sketch_med = nan(n_cases, n_specs);
angle_sketch_med = nan(n_cases, n_specs);

%% Main loop
for kk = 1:n_cases
    iiMax = iiMax_list(kk);
    U = Ufull(:, 1:iiMax);
    Q = orth(U);
    r = size(Q, 2);

    ell_values = [r, 2*r, 4*r, floor(n/2), n-r, n-1, n];
    % ell_values = [r,floor(n/2), n-(4*r), n-(2*r), n-r, n-1, n];
    valid = ell_values >= r & ell_values <= n;

    dimQ(kk) = r;
    sub_err(kk) = subspace(v_true, Q);
    spec_valid(kk, :) = valid;

    full_sol = solve_pencil(AA, Q, [], sketch_type, lam_true, v_true);
    abserr_full(kk) = full_sol.abserr;
    angle_full(kk) = full_sol.angle;

    sketch_err = nan(n_rep, n_specs);
    sketch_angle = nan(n_rep, n_specs);
    first_rep_sol = cell(1, n_specs);

    for rep = 1:n_rep
        for jj = 1:n_specs
            if ~valid(jj)
                continue;
            end

            ell = ell_values(jj);
            sol = solve_pencil(AA, Q, ell, sketch_type, lam_true, v_true);

            sketch_err(rep, jj) = sol.abserr;
            sketch_angle(rep, jj) = sol.angle;

            if rep == 1
                first_rep_sol{jj} = sol;
            end
        end

        if rep == 1 && ismember(r, plot_dims)
            plot_scan_compare(AA, Q, first_rep_sol, ell_labels, valid, lam_grid, r, rep);
        end
    end

    abserr_sketch_med(kk, :) = median(sketch_err, 1, 'omitnan');
    angle_sketch_med(kk, :) = median(sketch_angle, 1, 'omitnan');
end

%% diagnostics for l=20
Q = orth(Ufull(:, 1:iiMax_list(end)));
r = size(Q, 2);
A = AA;

lambda_true = lam_true;
x_true_red = Q' * v_true;
x_true_red = x_true_red / norm(x_true_red);

M0 = [Q, A*Q];
[~, ~, V0] = svd(M0, "econ", "vector");
Vtail0 = V0(:, r+1:2*r);
AQ = A * Q;

y_hat = Q' * v_true;
y_hat = y_hat / norm(y_hat);
y_star = [-lambda_true * y_hat; y_hat];
y_star = y_star / norm(y_star);
sin_y_star_V2 = norm(y_star - Vtail0 * (Vtail0' * y_star));

l_list = 2*r:n;
ntrials = 20;
gamma_reg = 1e-8;
gep_results = zeros(numel(l_list), 22);

for ii = 1:numel(l_list)
    l = l_list(ii);

    trial_metrics = zeros(ntrials, 7, 3);

    for tt = 1:ntrials
        G = randn(l, n);
        S_gaussian = G / sqrt(l);
        S_haar = orth(G')' * sqrt(n) / sqrt(l);
        srht_block = apply_srht([Q, AQ], l);

        SQ_list = {S_gaussian * Q, S_haar * Q, srht_block(:, 1:r)};
        SAQ_list = {S_gaussian * AQ, S_haar * AQ, srht_block(:, r+1:2*r)};

        for ss = 1:3
            Ms = [SQ_list{ss}, SAQ_list{ss}];
            [~, ~, Vs] = svd(Ms, "econ", "vector");

            Vhead_s = Vs(:, 1:r);
            % sin_tildeV2_V2_fro = norm(Vhead_s' * Vtail0, "fro");
            sin_tildeV2_V2_fro = norm(Vhead_s' * Vtail0);
            Vtail_s = Vs(:, r+1:2*r);
            sin_y_star_tildeV2 = norm(y_star - Vtail_s * (Vtail_s' * y_star));


            V11s = Vs(1:r, 1:r);
            V21s = Vs(r+1:2*r, 1:r);

            Ag = V21s';
            Bg = V11s';
            condV11 = cond(V11s);

            [Xg, Dg] = eig(Ag, Bg);
            lambda_g = diag(Dg);
            [eig_err, idx] = min(abs(lambda_g - lambda_true));

            lambda_sel = lambda_g(idx);
            x_sel = Xg(:, idx);
            x_sel = x_sel / norm(x_sel);

            q_sel = Q * x_sel;
            q_sel = q_sel / norm(q_sel);
            vec_angle = acos(min(1, abs(v_true' * q_sel)));

            [Yg, Dleft] = eig(Ag', Bg');
            lambda_left = diag(Dleft);
            [~, idx_left] = min(abs(lambda_left - conj(lambda_sel)));
            y_sel = Yg(:, idx_left);
            y_sel = y_sel / norm(y_sel);

            denom = abs(y_sel' * Bg * x_sel);
            numer = norm(x_sel) * norm(y_sel) * ...
                sqrt(norm(Ag, 2)^2 + abs(lambda_sel)^2 * norm(Bg, 2)^2);
            eig_cond = numer / denom;

            true_gep_residual = norm(Ag*x_true_red - lambda_true*Bg*x_true_red) / ...
                ((norm(Ag, 2) + abs(lambda_true)*norm(Bg, 2)) * norm(x_true_red));

            trial_metrics(tt, :, ss) = [
                sin_tildeV2_V2_fro, ...
                condV11, ...
                eig_cond, ...
                eig_err, ...
                vec_angle, ...
                true_gep_residual, ...
                sin_y_star_tildeV2
            ];
        end
    end

    gep_results(ii, :) = [
        l, ...
        mean(trial_metrics(:, :, 1), 1), ...
        mean(trial_metrics(:, :, 2), 1), ...
        mean(trial_metrics(:, :, 3), 1)
    ];
end

T = array2table(gep_results, ...
    "VariableNames", ...
    ["l", ...
     "gaussian_mean_sin_tildeV2_V2_fro", ...
     "gaussian_mean_cond_V11", ...
     "gaussian_mean_selected_GEP_eig_cond", ...
     "gaussian_mean_eig_error_GEP", ...
     "gaussian_mean_vector_angle", ...
     "gaussian_mean_true_solution_GEP_residual", ...
     "gaussian_mean_sin_y_star_tildeV2", ...
     "haar_mean_sin_tildeV2_V2_fro", ...
     "haar_mean_cond_V11", ...
     "haar_mean_selected_GEP_eig_cond", ...
     "haar_mean_eig_error_GEP", ...
     "haar_mean_vector_angle", ...
     "haar_mean_true_solution_GEP_residual", ...
     "haar_mean_sin_y_star_tildeV2", ...
     "srht_mean_sin_tildeV2_V2_fro", ...
     "srht_mean_cond_V11", ...
     "srht_mean_selected_GEP_eig_cond", ...
     "srht_mean_eig_error_GEP", ...
     "srht_mean_vector_angle", ...
     "srht_mean_true_solution_GEP_residual", ...
     "srht_mean_sin_y_star_tildeV2"]);

disp(T);
%%
figure;
plot(gep_results(:, 1), gep_results(:, 2), 'bo-', 'MarkerSize', 4,'LineWidth', 1.5);
hold on;
% plot(gep_results(:, 1), gep_results(:, 9), 'rs-', 'MarkerSize', 5,'LineWidth', 1.5);
plot(gep_results(:, 1), gep_results(:, 16), 'rs-', 'MarkerSize', 4,'LineWidth', 1.5);

l_axis = gep_results(:, 1);
y_data = gep_results(:, [2, 9, 16]);
ref_sqrt_n_minus_l = (1 ./ sqrt(l_axis)) .* sqrt(max(1 - l_axis ./ n, 0));
ref_inv_sqrt_l = 1 ./ sqrt(l_axis);
C_ref = max(y_data(:)) / max([ref_sqrt_n_minus_l; ref_inv_sqrt_l]);

% plot(l_axis, C_ref * ref_sqrt_n_minus_l, 'k--', 'LineWidth', 1.5);
% plot(l_axis, C_ref * ref_inv_sqrt_l, 'm--', 'LineWidth', 1.5);
xlabel('$\ell$', 'Interpreter', 'latex');
ylabel('$\|\widetilde V_1^* V_2\|_2$', 'Interpreter', 'latex');
% title('$\sin(\Theta)$ diagnostic for sketched right singular subspace', 'Interpreter', 'latex');
legend({'Gaussian', 'SRHT'}, ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', 13);
grid on;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 6.6, 5.8]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'saved_plots/testpara_gep_sin_tildeV2_V2.eps', '-depsc2', '-painters');

figure;
semilogy(gep_results(:, 1), gep_results(:, 5), 'bo-', 'MarkerSize', 4, 'LineWidth', 1.5);
hold on;
semilogy(gep_results(:, 1), gep_results(:, 6), 'b--o', 'MarkerSize', 4, 'LineWidth', 1.5);
semilogy(gep_results(:, 1), gep_results(:, 19), 'rs-', 'MarkerSize', 4, 'LineWidth', 1.5);
semilogy(gep_results(:, 1), gep_results(:, 20), 'r--s', 'MarkerSize', 4, 'LineWidth', 1.5);
xlabel('$\ell$', 'Interpreter', 'latex');
ylabel('error', 'Interpreter', 'latex');
legend({'Gaussian eig', 'Gaussian vector', 'SRHT eig', 'SRHT vector'}, ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', 13);
grid on;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 6.6, 5.8]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'saved_plots/testpara_gep_eig_vector_gaussian_srht.eps', '-depsc2', '-painters');

figure;
semilogy(gep_results(:, 1), gep_results(:, 8), 'bo-', 'MarkerSize', 4, 'LineWidth', 1.5);
hold on;
semilogy(gep_results(:, 1), gep_results(:, 22), 'rs-', 'MarkerSize', 4, 'LineWidth', 1.5);
semilogy(gep_results(:, 1), max(sin_y_star_V2, eps) * ones(size(gep_results(:, 1))), ...
    'k--', 'LineWidth', 1.5);
xlabel('$\ell$', 'Interpreter', 'latex');
ylabel('$\sin\angle(y_*, \mathrm{span}(\widetilde V_2))$', 'Interpreter', 'latex');
legend({'Gaussian', 'SRHT', 'Original $V_2$'}, ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', 13);
grid on;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 6.6, 5.8]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'saved_plots/testpara_gep_y_star_tildeV2_angle.eps', '-depsc2', '-painters');
%%
plot_summary(dimQ, sub_err, abserr_full, angle_full, ...
    abserr_sketch_med, angle_sketch_med, spec_valid, ell_labels, sketch_type);

%% Helpers
function sol = solve_pencil(A, Q, ell, sketch_type, lam_true, v_true)
    AQ = A * Q;

    if isempty(ell)
        Asketch = AQ;
        Qsketch = Q;
    else
        [Qsketch, Asketch] = build_sketch_blocks(Q, AQ, ell, sketch_type);
    end

    if size(Asketch, 1) > size(Asketch, 2)
        [X, lam_all] = local_tls_pencil_eigs(Asketch, Qsketch);
    elseif size(Asketch, 1) == size(Asketch, 2)
        [X, lam_all] = eig(Asketch, Qsketch, 'vector');
    else
        error('Fat pencils are not supported: ell = %d, r = %d.', size(Asketch, 1), size(Asketch, 2));
    end

    r = numel(lam_all);
    lam_rq = nan(r, 1);
    Uritz = Q * X;

    for ii = 1:r
        u = Uritz(:, ii);
        u = u / norm(u);
        lam_rq(ii) = (u' * A * u) / (u' * u);
        Uritz(:, ii) = u;
    end

    score = abs(lam_rq - lam_true);
    score(~isfinite(score)) = inf;
    [~, idx] = min(score);

    u = Uritz(:, idx);

    sol = struct();
    sol.Asketch = Asketch;
    sol.Qsketch = Qsketch;
    sol.lam = lam_rq(idx);
    sol.abserr = abs(sol.lam - lam_true);
    sol.angle = subspace(v_true, u);
end

function [Qsketch, Asketch] = build_sketch_blocks(Q, AQ, ell, sketch_type)
    [n, r] = size(Q);
    sketch_type = lower(char(sketch_type));

    if strcmp(sketch_type, 'gaussian')
        G = randn(ell, n);
        S = G / sqrt(ell);
        Qsketch = S * Q;
        Asketch = S * AQ;

    elseif strcmp(sketch_type, 'haar')
        G = randn(ell, n);
        S = orth(G')' * sqrt(n) / sqrt(ell);
        Qsketch = S * Q;
        Asketch = S * AQ;

    elseif strcmp(sketch_type, 'srht')
        sketched = apply_srht([Q, AQ], ell);
        Qsketch = sketched(:, 1:r);
        Asketch = sketched(:, r+1:2*r);

    else
        error('Unknown sketch_type "%s". Use gaussian, haar, or srht.', sketch_type);
    end
end

function [vec_tls_all, lam_tls_all] = local_tls_pencil_eigs(A, B)
    [nA, rA] = size(A);
    [nB, rB] = size(B);
    if nA ~= nB || rA ~= rB
        error('A and B must have the same size n-by-r.');
    end

    r = rA;
    [~, ~, Vc] = svd([B, A], 0);
    if size(Vc, 2) < r
        error('SVD did not return enough right singular vectors.');
    end

    V11 = Vc(1:r, 1:r);
    V21 = Vc(r+1:2*r, 1:r);
    [vec_tls_all, lam_tls_all] = eig(V21', V11', 'vector');
end

function s = sigma_min_pencil(A, B, lam)
    svals = svd(A - lam * B, 'econ');
    s = svals(end);
end

function SX = apply_srht(X, ell)
    n = size(X, 1);
    p = 2^nextpow2(n);

    Xpad = zeros(p, size(X, 2));
    Xpad(1:n, :) = X;

    signs = 2 * (rand(p, 1) > 0.5) - 1;
    Xpad = signs .* Xpad;

    HX = fast_hadamard_transform(Xpad) / sqrt(p);
    rows = randperm(p, ell);
    SX = sqrt(p / ell) * HX(rows, :);
end

function Y = fast_hadamard_transform(Y)
    p = size(Y, 1);
    h = 1;

    while h < p
        for first = 1:2*h:p
            top_idx = first:first+h-1;
            bot_idx = first+h:first+2*h-1;

            top = Y(top_idx, :);
            bot = Y(bot_idx, :);

            Y(top_idx, :) = top + bot;
            Y(bot_idx, :) = top - bot;
        end
        h = 2*h;
    end
end

function plot_scan_compare(A, Q, sketch_solutions, labels, valid, lam_grid, r, rep)
    AQ = A * Q;
    sigma_full = nan(size(lam_grid));
    valid_idx = find(valid);
    sigma_sketch = nan(numel(lam_grid), numel(valid_idx));

    for ii = 1:numel(lam_grid)
        lam = lam_grid(ii);
        sigma_full(ii) = sigma_min_pencil(AQ, Q, lam);
        for jj = 1:numel(valid_idx)
            sol = sketch_solutions{valid_idx(jj)};
            sigma_sketch(ii, jj) = sigma_min_pencil(sol.Asketch, sol.Qsketch, lam);
        end
    end

    figure;
    h = gobjects(numel(valid_idx) + 1, 1);
    h(1) = semilogy(lam_grid, sigma_full, 'b-', 'LineWidth', 1.5);
    hold on;

    colors = lines(max(numel(valid_idx), 1));
    for jj = 1:numel(valid_idx)
        h(jj+1) = semilogy(lam_grid, sigma_sketch(:, jj), '-', ...
            'Color', colors(jj, :), 'LineWidth', 1.5);
    end

    xlabel('\lambda');
    ylabel('\sigma_{min}');
    legend(h, [{'unsketched'}, labels(valid_idx)], 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6.6, 5.8]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, sprintf('saved_plots/testpara_curve04_%d.eps', r), '-depsc2', '-painters');
end

function plot_summary(dimQ, sub_err, abserr_full, angle_full, ...
    abserr_sketch_med, angle_sketch_med, spec_valid, labels, sketch_type)

    valid_cols = find(any(spec_valid, 1));
    colors = lines(max(numel(valid_cols), 1));

    figure;
    h = gobjects(numel(valid_cols) + 2, 1);
    h(1) = semilogy(dimQ, angle_full, 'bo-', 'LineWidth', 1.5);
    hold on;

    for jj = 1:numel(valid_cols)
        col = valid_cols(jj);
        h(jj+1) = semilogy(dimQ, angle_sketch_med(:, col), '-', ...
            'Color', colors(jj, :), 'LineWidth', 1.5);
    end

    h(end) = semilogy(dimQ, sub_err, 'k--', 'LineWidth', 1.5);
    xlabel('dim(Q)', 'Interpreter', 'latex');
    legend(h, [{'MP'}, labels(valid_cols), {'subspace angle'}], 'Location', 'best');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6.6, 5.8]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, 'saved_plots/testpara_angle04.eps', '-depsc2', '-painters');

    figure;
    h = gobjects(numel(valid_cols) + 2, 1);
    h_idx = 1;
    h(h_idx) = semilogy(dimQ, abserr_full, 'bo-', 'LineWidth', 1.5);
    hold on;

    err_labels = cell(1, numel(h));
    err_labels{1} = 'MP refined';
    label_idx = 2;

    for jj = 1:numel(valid_cols)
        col = valid_cols(jj);
        h_idx = h_idx + 1;
        h(h_idx) = semilogy(dimQ, abserr_sketch_med(:, col), '-', ...
            'Color', colors(jj, :), 'LineWidth', 1.5);
        err_labels{label_idx} = [labels{col}, ' refined'];
        label_idx = label_idx + 1;
    end

    h_idx = h_idx + 1;
    h(h_idx) = semilogy(dimQ, sub_err, 'k--', 'LineWidth', 1.5);
    err_labels{h_idx} = 'subspace angle';

    xlabel('dim(Q)', 'Interpreter', 'latex');
    legend(h(1:h_idx), err_labels(1:h_idx), 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6.6, 5.8]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, 'saved_plots/testpara_err04.eps', '-depsc2', '-painters');

    figure;
    semilogy(dimQ, abserr_full, 'b-o', 'LineWidth', 1.5);
    hold on;
    semilogy(dimQ, angle_full, 'b--o', 'LineWidth', 1.5);

    for jj = 1:numel(valid_cols)
        col = valid_cols(jj);
        semilogy(dimQ, abserr_sketch_med(:, col), '-', ...
            'Color', colors(jj, :), 'LineWidth', 1.5);

        semilogy(dimQ, angle_sketch_med(:, col), '--', ...
            'Color', colors(jj, :), 'LineWidth', 1.5);
    end
    semilogy(dimQ, sub_err, 'k--', 'LineWidth', 1.5);

    legend_handles = gobjects(numel(valid_cols) + 2, 1);
    legend_labels = [{'MP'}, labels(valid_cols), {'subspace angle'}];

    legend_handles(1) = plot(nan, nan, 'b-', 'LineWidth', 1.5);
    for jj = 1:numel(valid_cols)
        legend_handles(jj+1) = plot(nan, nan, '-', ...
            'Color', colors(jj, :), 'LineWidth', 1.5);
    end
    legend_handles(end) = plot(nan, nan, 'k--', 'LineWidth', 1.5);

    xlabel('dim(Q)', 'Interpreter', 'latex');
    ylabel('error', 'Interpreter', 'latex');
    lgd = legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 12);
    set(lgd, 'Box', 'off', 'Color', 'none');
    grid on;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6.6, 5.8]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPositionMode', 'auto');

    print(gcf, 'saved_plots/testpara_err_angle_combined04.eps', '-depsc2', '-painters');
end

function report = resolve_local_cluster(problem, sigma_fun_scan, sigma_fun_refine, ...
        cand, global_interval, refine_opts, cfg, cluster_opts)
%RESOLVE_LOCAL_CLUSTER  Local cluster-resolution post-processing for a candidate.
%   At a refined candidate lambda*, check whether the singular values of
%   Q_B(lambda*) are clustered. If so, estimate a local scale from a
%   finite-difference approximation of Q_B'(lambda*), build a finer local grid
%   around the original bracket, and rerun local-minimum detection plus
%   refinement to resolve hidden nearby dips (thin H-shape case).
%
%   QR policy follows cfg: the cluster singular-value check and the
%   finite-difference reference use the CPQR/truncated scan objective, while the
%   micro-scan refinement uses sigma_fun_refine (plain, untruncated QR).

    lam_star = cand.refined_lambda;
    L_macro = max(abs(cand.bracket - lam_star));

    [QB_star, factor_info] = build_reference_QB(problem, lam_star, cfg);
    [~, S] = svd(QB_star, 'econ');
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
    report.response = [];
    report.refined_sigma_original = cand.refined_sigma;
    report.merged_lambda = lam_star;
    report.merged_sigma = cand.refined_sigma;
    report.merged_replaced = false;
    report.n_merged = 0;
    report.dedup_tol = [];
    report.reference_rank = factor_info.rank;
    report.reference_cols = factor_info.selected_cols;
    report.micro_scan = struct('lamvec', [], 'S', [], 'candidate_idx', []);
    report.micro_summary = struct('refined_candidates', [], 'extra_mask', []);
    report.message = '';

    if ~report.cluster_detected
        report.message = 'No cluster detected at the refined lambda.';
        return;
    end

    pair_opts = struct();
    pair_opts.normalize_columns = false;
    pair_opts.qr_tau = get_cfg_value(cfg, 'qr_tau', []);
    pair_opts.pivot = get_cfg_flag(cfg, 'qr_pivot');
    pair_opts.sign_fix = true;

    [dQB, ~, ~, deriv_info] = approx_QB_derivative_from_Aop( ...
        problem.ops.A, problem.meta.mB, lam_star, cluster_opts.fd_step, pair_opts);
    response = norm(dQB, 2);
    % response = 5e-2;  % manual shortcut; keep the adaptive response above active
    report.response = response;
    report.derivative_info = deriv_info;

    if response <= eps(class(response))
        report.message = 'Cluster detected, but the finite-difference response is numerically zero.';
        return;
    end

    L_fine = 0.5 * sigma_min / response;
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
        S_micro(i) = sigma_fun_scan(lamvec(i));
    end

    J = 2:numel(lamvec)-1;
    J = J(S_micro(J) < S_micro(J-1) & S_micro(J) < S_micro(J+1));

    micro_scan = struct();
    micro_scan.lamvec = lamvec;
    micro_scan.S = S_micro;
    micro_scan.candidate_idx = J;

    refined_candidates = refine_candidates(sigma_fun_refine, micro_scan, refine_opts);

    % Dedup tolerance must exceed the minimizer's positional jitter on the flat
    % clustered valley, which spans several micro-grid steps; L_fine alone is too
    % tight and re-flags the original dip as a phantom extra.
    dedup_factor = get_cfg_value(cluster_opts, 'dedup_factor', 8);
    lam_tol = dedup_factor * L_fine;
    [extra_mask, dup_mask] = local_extra_dip_mask(refined_candidates, lam_star, lam_tol);

    % Merge the original dip with its micro re-detections (within lam_tol): if any
    % re-detection is deeper, it replaces the original as this dip's localization.
    dup_cands = refined_candidates(dup_mask);
    report.n_merged = numel(dup_cands);
    for i = 1:numel(dup_cands)
        if dup_cands(i).refined_sigma < report.merged_sigma
            report.merged_lambda = dup_cands(i).refined_lambda;
            report.merged_sigma = dup_cands(i).refined_sigma;
            report.merged_replaced = true;
        end
    end

    report.dedup_tol = lam_tol;
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


function [extra_mask, dup_mask] = local_extra_dip_mask(candidates, lam_star, lam_tol)
    n = numel(candidates);
    extra_mask = false(n, 1);
    dup_mask = false(n, 1);
    for i = 1:n
        if ~candidates(i).accepted
            continue;
        end
        if abs(candidates(i).refined_lambda - lam_star) > lam_tol
            extra_mask(i) = true;   % genuinely new dip
        else
            dup_mask(i) = true;     % re-detection of the original dip
        end
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

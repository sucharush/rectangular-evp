function candidates = refine_candidates(sigma_fun, scan, opts)
% Refine each scan candidate into a candidate record.
%
% Required fields in opts:
%   opts.bracket_halfwidth
%   opts.sigma_cut
%   opts.minimizer
%   opts.minimizer_opts
%
% Optional:
%   opts.pre_refine_filter   = function handle @(cand, scan, opts) -> true/false
%   opts.local_analyzer      = function handle @(sigma_fun, cand, opts) -> cand
%
% Output:
%   candidates : struct array, one per detected scan dip

    required_fields = {'bracket_halfwidth', 'sigma_cut', 'minimizer', 'minimizer_opts'};
    for k = 1:numel(required_fields)
        if ~isfield(opts, required_fields{k})
            error('refine_candidates: missing opts.%s', required_fields{k});
        end
    end

    if ~isfield(opts, 'pre_refine_filter')
        opts.pre_refine_filter = [];
    end
    if ~isfield(opts, 'local_analyzer')
        opts.local_analyzer = [];
    end

    lamvec = scan.lamvec;
    S = scan.S;
    J = scan.candidate_idx;

    candidates = repmat(empty_candidate_template(), numel(J), 1);

    for t = 1:numel(J)
        j = J(t);

        jl = max(1, j - opts.bracket_halfwidth);
        jr = min(numel(lamvec), j + opts.bracket_halfwidth);

        a = lamvec(jl);
        b = lamvec(jr);

        cand = empty_candidate_template();

        % -------- raw scan information --------
        cand.id = t;
        cand.scan_index = j;
        cand.scan_lambda = lamvec(j);
        cand.scan_sigma = S(j);

        cand.left_index = jl;
        cand.right_index = jr;
        cand.bracket = [a, b];

        cand.status = 'raw';
        cand.accepted = false;

        cand.logs = struct();
        cand.logs.n_calls = 0;
        cand.logs.message = '';

        cand.diagnostics = struct();

        % -------- optional pre-filter --------
        if ~isempty(opts.pre_refine_filter)
            keep = opts.pre_refine_filter(cand, scan, opts);
            if ~keep
                cand.status = 'rejected_prefilter';
                cand.logs.message = 'Rejected by pre_refine_filter.';
                candidates(t) = cand;
                continue;
            end
        end

        % -------- local refinement --------
        % counted_obj = make_counted_objective(@(lam) problem.ops.sigma(lam));
        counted_obj = make_counted_objective(@(lam) sigma_fun(lam));

        try
            [lam_star, sig_star, info] = opts.minimizer(a, b, counted_obj, opts.minimizer_opts);

            cand.refined_lambda = lam_star;
            cand.refined_sigma = sig_star;
            cand.minimizer_info = info;
            cand.logs.n_calls = counted_obj.get_count();

            if ~isfinite(lam_star) || ~isfinite(sig_star)
                cand.status = 'failed_nonfinite';
                cand.logs.message = 'Minimizer returned non-finite values.';
            elseif sig_star > opts.sigma_cut
                cand.status = 'rejected_sigma_cut';
                cand.logs.message = sprintf('Rejected by sigma_cut: %.3e > %.3e', sig_star, opts.sigma_cut);
            else
                cand.status = 'accepted';
                cand.accepted = true;
            end

        catch ME
            cand.status = 'failed_exception';
            cand.logs.n_calls = counted_obj.get_count();
            cand.logs.message = ME.message;
            cand.minimizer_info = struct('status', 'exception');
        end

        % -------- optional local analyzer --------
        if ~isempty(opts.local_analyzer)
            try
                % cand = opts.local_analyzer(problem, cand, opts);
                cand = opts.local_analyzer(sigma_fun, cand, opts);
            catch ME
                if isempty(cand.logs.message)
                    cand.logs.message = sprintf('local_analyzer failed: %s', ME.message);
                else
                    cand.logs.message = sprintf('%s | local_analyzer failed: %s', cand.logs.message, ME.message);
                end
            end
        end

        candidates(t) = cand;
    end
end


function cand = empty_candidate_template()
    cand = struct();

    % identity / provenance
    cand.id = [];
    cand.scan_index = [];
    cand.scan_lambda = [];
    cand.scan_sigma = [];

    cand.left_index = [];
    cand.right_index = [];
    cand.bracket = [];

    % refined result
    cand.refined_lambda = [];
    cand.refined_sigma = [];
    cand.minimizer_info = struct();

    % status
    cand.status = '';
    cand.accepted = false;

    % bookkeeping
    cand.logs = struct();
    cand.diagnostics = struct();
end

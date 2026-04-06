% function result = solve_scan_refine(problem, opts)
% % General scan-refine solver based on problem.ops.sigma(lambda).
% %
% % Required:
% %   opts.scan
% %   opts.refine
% %
% % Output:
% %   result.scan
% %   result.candidates
% %   result.summary
% %   result.logs
% 
%     if ~isfield(opts, 'scan')
%         error('solve_scan_refine: opts.scan is required.');
%     end
%     if ~isfield(opts, 'refine')
%         error('solve_scan_refine: opts.refine is required.');
%     end
% 
%     result = struct();
%     result.method = 'scan_refine';
%     result.problem_name = '';
%     result.problem_family = '';
% 
%     if isfield(problem, 'name')
%         result.problem_name = problem.name;
%     end
%     if isfield(problem, 'family')
%         result.problem_family = problem.family;
%     end
% 
%     % -------- stage 1: global scan --------
%     scan = scan_sigma(problem, opts.scan);
% 
%     % -------- stage 2: candidate refinement --------
%     candidates = refine_candidates(problem, scan, opts.refine);
% 
%     % -------- stage 3: compact summary --------
%     summary = summarize_candidates(candidates);
% 
%     % -------- logs --------
%     logs = struct();
%     logs.n_scan_points = numel(scan.lamvec);
%     logs.n_raw_candidates = numel(scan.candidate_idx);
%     logs.n_total_candidates = numel(candidates);
%     logs.n_accepted_candidates = summary.n_accepted;
% 
%     if isempty(candidates)
%         logs.total_refine_calls = 0;
%     else
%         n_calls = arrayfun(@(c) get_log_calls(c), candidates);
%         logs.total_refine_calls = sum(n_calls);
%     end
% 
%     result.scan = scan;
%     result.candidates = candidates;
%     result.summary = summary;
%     result.logs = logs;
% end
% 
% 
% function n = get_log_calls(cand)
%     n = 0;
%     if isfield(cand, 'logs') && isfield(cand.logs, 'n_calls') && ~isempty(cand.logs.n_calls)
%         n = cand.logs.n_calls;
%     end
% end

function result = solve_scan_refine(sigma_fun, opts)
% General scan-refine solver for a scalar objective sigma_fun(lambda).
%
% Required:
%   sigma_fun
%   opts.scan
%   opts.refine

    if nargin < 2
        error('solve_scan_refine: two inputs required: sigma_fun, opts.');
    end
    if ~isa(sigma_fun, 'function_handle')
        error('solve_scan_refine: sigma_fun must be a function handle.');
    end
    if ~isfield(opts, 'scan')
        error('solve_scan_refine: opts.scan is required.');
    end
    if ~isfield(opts, 'refine')
        error('solve_scan_refine: opts.refine is required.');
    end

    scan = scan_sigma(sigma_fun, opts.scan);
    candidates = refine_candidates(sigma_fun, scan, opts.refine);
    summary = summarize_candidates(candidates);

    logs = struct();
    logs.n_scan_points = numel(scan.lamvec);
    logs.n_raw_candidates = numel(scan.candidate_idx);
    logs.n_total_candidates = numel(candidates);
    logs.n_accepted_candidates = summary.n_accepted;

    if isempty(candidates)
        logs.total_refine_calls = 0;
    else
        n_calls = zeros(numel(candidates),1);
        for k = 1:numel(candidates)
            if isfield(candidates(k), 'logs') && isfield(candidates(k).logs, 'n_calls') ...
                    && ~isempty(candidates(k).logs.n_calls)
                n_calls(k) = candidates(k).logs.n_calls;
            end
        end
        logs.total_refine_calls = sum(n_calls);
    end

    result = struct();
    result.method = 'scan_refine';
    result.scan = scan;
    result.candidates = candidates;
    result.summary = summary;
    result.logs = logs;
end
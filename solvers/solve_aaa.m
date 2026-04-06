% function result = solve_aaa(problem, opts)
% % General AAA-based rational surrogate construction for a matrix-valued
% % operator F(lambda) on a real interval [a,b].
% %
% % Required:
% %   problem.ops.F
% %   opts.interval = [a, b]
% %   opts.aaa      = struct()
% %
% % Optional:
% %   opts.verify.do_scan = false/true
% %   opts.verify.lam_grid
% %   opts.verify.use_problem_sigma = false/true
% %
% % Output:
% %   result.raw
% %   result.summary
% %   result.verification
% %   result.logs
% 
%     % ----------------------------
%     % checks
%     % ----------------------------
%     if ~isfield(problem, 'ops') || ~isfield(problem.ops, 'F')
%         error('solve_aaa: problem.ops.F is required.');
%     end
% 
%     if ~isfield(opts, 'interval') || numel(opts.interval) ~= 2
%         error('solve_aaa: opts.interval = [a,b] is required.');
%     end
% 
%     if ~isfield(opts, 'aaa')
%         error('solve_aaa: opts.aaa is required.');
%     end
% 
%     if ~isfield(opts, 'verify')
%         opts.verify = struct();
%     end
%     if ~isfield(opts.verify, 'do_scan')
%         opts.verify.do_scan = false;
%     end
%     if ~isfield(opts.verify, 'use_problem_sigma')
%         opts.verify.use_problem_sigma = false;
%     end
% 
%     a = opts.interval(1);
%     b = opts.interval(2);
%     F_op = problem.ops.F;
% 
%     % ----------------------------
%     % construct surrogate / pencil
%     % ----------------------------
%     raw = aaa_rect_pencil_interval(F_op, a, b, opts.aaa);
% 
%     % ----------------------------
%     % summary
%     % ----------------------------
%     summary = struct();
%     summary.interval = [a, b];
%     summary.degree = raw.m;
%     summary.n_support_points = numel(raw.zj);
%     summary.support_points = raw.zj;
%     summary.weights = raw.wj;
%     summary.sample_count = numel(raw.Z);
% 
%     % ----------------------------
%     % optional verification
%     % ----------------------------
%     verification = struct();
% 
%     if opts.verify.do_scan
%         if ~isfield(opts.verify, 'lam_grid') || isempty(opts.verify.lam_grid)
%             error('solve_aaa: opts.verify.lam_grid is required when verify.do_scan=true.');
%         end
% 
%         lam_grid = opts.verify.lam_grid(:);
%         sigma_original = zeros(size(lam_grid));
%         sigma_bary = zeros(size(lam_grid));
%         sigma_newton = zeros(size(lam_grid));
% 
%         for i = 1:numel(lam_grid)
%             lam = lam_grid(i);
% 
%             M0 = F_op(lam);
%             Mb = aaa_eval_matrix_barycentric(lam, raw.zj, raw.wj, raw.D);
%             Mn = aaa_eval_matrix_newton(lam, raw.D, raw.sigma, raw.beta, raw.h, raw.k);
% 
%             if opts.verify.use_problem_sigma && isfield(problem.ops, 'sigma') && ~isempty(problem.ops.sigma)
%                 sigma_original(i) = problem.ops.sigma(lam);
%             else
%                 sigma_original(i) = local_target_sigma(M0, raw);
%             end
% 
%             sigma_bary(i) = local_target_sigma(Mb, raw);
%             sigma_newton(i) = local_target_sigma(Mn, raw);
%         end
% 
%         verification.lam_grid = lam_grid;
%         verification.sigma_original = sigma_original;
%         verification.sigma_bary = sigma_bary;
%         verification.sigma_newton = sigma_newton;
%     end
% 
%     % ----------------------------
%     % logs
%     % ----------------------------
%     logs = struct();
%     logs.interval = [a, b];
%     logs.degree = raw.m;
%     logs.n_support_points = numel(raw.zj);
% 
%     % ----------------------------
%     % output
%     % ----------------------------
%     result = struct();
%     result.method = 'aaa';
%     result.problem_name = '';
%     result.problem_family = '';
% 
%     if isfield(problem, 'name'),   result.problem_name = problem.name; end
%     if isfield(problem, 'family'), result.problem_family = problem.family; end
% 
%     result.raw = raw;
%     result.summary = summary;
%     result.verification = verification;
%     result.logs = logs;
% end
% 
% 
% function s = local_target_sigma(M, raw)
% % Default singular-value tracking rule for verification.
% %
% % For a p-by-q operator with p >= q, the natural default is sigma_q,
% % i.e. the smallest singular value of the nonzero band if full column rank.
% %
% % If raw.q is unavailable, fall back to the smallest singular value.
% 
%     svals = svd(M, 'econ');
% 
%     if isempty(svals)
%         s = NaN;
%         return;
%     end
% 
%     if isfield(raw, 'q') && ~isempty(raw.q) && raw.q <= numel(svals)
%         s = svals(raw.q);
%     else
%         s = svals(end);
%     end
% end
% 
function result = solve_aaa(F_op, opts)
% General AAA-based rational surrogate construction for a matrix-valued
% operator F_op(lambda) on a real interval [a,b].
%
% Required:
%   F_op
%   opts.interval = [a, b]
%   opts.aaa      = struct()
%
% Optional:
%   opts.verify.do_scan = false/true
%   opts.verify.lam_grid
%   opts.verify.sigma_fun = []   % optional scalar reference built from F_op
%
% Output:
%   result.raw
%   result.summary
%   result.verification
%   result.logs

    % ----------------------------
    % checks
    % ----------------------------
    if nargin < 2
        error('solve_aaa: two inputs required: F_op, opts.');
    end
    if ~isa(F_op, 'function_handle')
        error('solve_aaa: F_op must be a function handle.');
    end

    if ~isfield(opts, 'interval') || numel(opts.interval) ~= 2
        error('solve_aaa: opts.interval = [a,b] is required.');
    end

    if ~isfield(opts, 'aaa')
        error('solve_aaa: opts.aaa is required.');
    end

    if ~isfield(opts, 'verify')
        opts.verify = struct();
    end
    if ~isfield(opts.verify, 'do_scan')
        opts.verify.do_scan = false;
    end
    if ~isfield(opts.verify, 'sigma_fun')
        opts.verify.sigma_fun = [];
    end

    a = opts.interval(1);
    b = opts.interval(2);

    % ----------------------------
    % construct surrogate / pencil
    % ----------------------------
    raw = aaa_rect_pencil_interval(F_op, a, b, opts.aaa);

    % ----------------------------
    % summary
    % ----------------------------
    summary = struct();
    summary.interval = [a, b];
    summary.degree = raw.m;
    summary.n_support_points = numel(raw.zj);
    summary.support_points = raw.zj;
    summary.weights = raw.wj;
    summary.sample_count = numel(raw.Z);

    % ----------------------------
    % optional verification
    % ----------------------------
    verification = struct();

    if opts.verify.do_scan
        if ~isfield(opts.verify, 'lam_grid') || isempty(opts.verify.lam_grid)
            error('solve_aaa: opts.verify.lam_grid is required when verify.do_scan=true.');
        end

        lam_grid = opts.verify.lam_grid(:);
        sigma_original = zeros(size(lam_grid));
        sigma_bary = zeros(size(lam_grid));
        sigma_newton = zeros(size(lam_grid));

        for i = 1:numel(lam_grid)
            lam = lam_grid(i);

            % original matrix-valued operator
            M0 = F_op(lam);

            % barycentric / Newton surrogates
            Mb = aaa_eval_matrix_barycentric(lam, raw.zj, raw.wj, raw.D);
            Mn = aaa_eval_matrix_newton(lam, raw.D, raw.sigma, raw.beta, raw.h, raw.k);

            % verification reference always comes from F_op directly
            if ~isempty(opts.verify.sigma_fun)
                sigma_original(i) = opts.verify.sigma_fun(lam);
            else
                sigma_original(i) = local_target_sigma(M0, raw);
            end

            sigma_bary(i) = local_target_sigma(Mb, raw);
            sigma_newton(i) = local_target_sigma(Mn, raw);
        end

        verification.lam_grid = lam_grid;
        verification.sigma_original = sigma_original;
        verification.sigma_bary = sigma_bary;
        verification.sigma_newton = sigma_newton;
    end

    % ----------------------------
    % logs
    % ----------------------------
    logs = struct();
    logs.interval = [a, b];
    logs.degree = raw.m;
    logs.n_support_points = numel(raw.zj);

    % ----------------------------
    % output
    % ----------------------------
    result = struct();
    result.method = 'aaa';
    result.raw = raw;
    result.summary = summary;
    result.verification = verification;
    result.logs = logs;
end


function s = local_target_sigma(M, raw)
% Default tracked singular value for verification.
% If raw.q is available, track sigma_q; otherwise use the smallest singular value.

    svals = svd(M, 'econ');

    if isempty(svals)
        s = NaN;
        return;
    end

    if isfield(raw, 'q') && ~isempty(raw.q) && raw.q <= numel(svals)
        s = svals(raw.q);
    else
        s = svals(end);
    end
end
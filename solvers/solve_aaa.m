% function result = solve_aaa(F_op, opts)
% % General AAA-based rational surrogate construction for a matrix-valued
% % operator F_op(lambda) on a real interval [a,b].
% %
% % Required:
% %   F_op
% %   opts.interval = [a, b]
% %   opts.aaa      = struct()
% %
% % Optional:
% %   opts.verify.do_scan = false/true
% %   opts.verify.lam_grid
% %   opts.verify.sigma_fun = []   % optional scalar reference built from F_op
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
%     if nargin < 2
%         error('solve_aaa: two inputs required: F_op, opts.');
%     end
%     if ~isa(F_op, 'function_handle')
%         error('solve_aaa: F_op must be a function handle.');
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
%     if ~isfield(opts.verify, 'sigma_fun')
%         opts.verify.sigma_fun = [];
%     end
% 
%     a = opts.interval(1);
%     b = opts.interval(2);
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
%         err_bary_fro = zeros(size(lam_grid));
%         err_bary_2 = zeros(size(lam_grid));
% 
%         for i = 1:numel(lam_grid)
%             lam = lam_grid(i);
% 
%             % original matrix-valued operator
%             M0 = F_op(lam);
% 
%             % barycentric / Newton surrogates
%             Mb = aaa_eval_matrix_barycentric(lam, raw.zj, raw.wj, raw.D);
%             Mn = aaa_eval_matrix_newton(lam, raw.D, raw.sigma, raw.beta, raw.h, raw.k);
% 
%             % verification reference always comes from F_op directly
%             if ~isempty(opts.verify.sigma_fun)
%                 sigma_original(i) = opts.verify.sigma_fun(lam);
%             else
%                 sigma_original(i) = local_target_sigma(M0, raw);
%             end
% 
%             sigma_bary(i) = local_target_sigma(Mb, raw);
%             sigma_newton(i) = local_target_sigma(Mn, raw);
%             err_bary_fro(i) = norm(M0 - Mb, 'fro')/norm(M0, 'fro');
%             err_bary_2(i) = norm(M0 - Mb)/norm(M0);
%         end
% 
%         verification.lam_grid = lam_grid;
%         verification.sigma_original = sigma_original;
%         verification.sigma_bary = sigma_bary;
%         verification.sigma_newton = sigma_newton;
%         verification.err_bary_fro = err_bary_fro;
%         verification.err_bary_2 = err_bary_2;
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
%     result.raw = raw;
%     result.summary = summary;
%     result.verification = verification;
%     result.logs = logs;
% end
% 
% 
% function s = local_target_sigma(M, raw)
% % Default tracked singular value for verification.
% % If raw.q is available, track sigma_q; otherwise use the smallest singular value.
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

function result = solve_aaa(op, opts)
% Build AAA surrogate and record barycentric reconstruction error
% at every sampled point Z(j).
%
% INPUT
%   op             : function handle
%   opts.interval  : [a, b]
%   opts.aaa       : AAA options
%
% Modes:
%   opts.aaa.method = 'direct'
%       op is QB_op(lambda)
%
%   opts.aaa.method = 'procrustes_from_A'
%       op is A_op(lambda), AAA routine internally builds aligned QB_all
%
% OUTPUT
%   result.raw
%   result.summary
%   result.verification
%   result.logs

    if nargin < 2
        error('solve_aaa: two inputs required: op, opts.');
    end
    if ~isa(op, 'function_handle')
        error('solve_aaa: op must be a function handle.');
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
    if ~isfield(opts.verify, 'compute_bary_error')
        opts.verify.compute_bary_error = true;
    end

    a = opts.interval(1);
    b = opts.interval(2);
    
    % t0 = tic;
    raw = aaa_rect_pencil_interval(op, a, b, opts.aaa);
    % elapsed_time = toc(t0);

    summary = struct();
    summary.interval = [a, b];
    summary.method = raw.method;
    summary.degree = raw.m;
    summary.n_support_points = numel(raw.zj);
    summary.support_points = raw.zj;
    summary.weights = raw.wj;
    summary.sample_count = numel(raw.Z);

    verification = struct();

    if opts.verify.compute_bary_error
        Z = raw.Z(:);
        nZ = numel(Z);

        relerr_bary = zeros(nZ, 1);
        abserr_bary = zeros(nZ, 1);

        for j = 1:nZ
            lam = Z(j);

            % exact sampled matrix at this grid point
            if strcmpi(raw.method, 'direct')
                Qtrue = op(lam);
            elseif strcmpi(raw.method, 'procrustes_from_a')
                Qtrue = raw.QB_all(:,:,j);
            else
                error('solve_aaa: unknown raw.method "%s".', raw.method);
            end

            Rb = aaa_eval_matrix_barycentric(lam, raw.zj, raw.wj, raw.D);
            % Rb = aaa_eval_matrix_newton(lam, raw.D, raw.sigma, raw.beta, raw.h, raw.k);
            % disp(Rb)

            abserr_bary(j) = norm(Qtrue - Rb, 'fro');
            % disp(abserr_bary(j))
            relerr_bary(j) = abserr_bary(j) / max(norm(Qtrue, 'fro'), 1e-14);
        end

        verification.Z = Z;
        verification.abserr_bary = abserr_bary;
        verification.relerr_bary = relerr_bary;
        verification.mean_relerr_bary = mean(relerr_bary);
        verification.max_relerr_bary = max(relerr_bary);
    end

    logs = struct();
    logs.interval = [a, b];
    logs.method = raw.method;
    logs.degree = raw.m;
    logs.elapsed_time = raw.elapsed_time;
    logs.n_support_points = numel(raw.zj);

    result = struct();
    result.method = 'aaa';
    result.raw = raw;
    result.summary = summary;
    result.verification = verification;
    result.logs = logs;
end


function s = local_target_sigma(M, raw)
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


function M = local_build_aligned_block_at_lambda(A_op, lam, raw)
% Lightweight reconstruction for verification at a single lambda.
% Uses the stored previous-alignment reference only approximately if needed.
% For exact reproducibility, prefer verifying on raw.Z sample points.

    A = A_op(lam);
    [Q, R] = qr(A, 0);

    if isfield(raw, 'sign_fix') && raw.sign_fix
        d = sign(diag(R));
        d(d == 0) = 1;
        Q = Q * diag(d);
    end

    if isfield(raw, 'mB')
        M = Q(1:raw.mB, :);
    else
        error('local_build_aligned_block_at_lambda: raw.mB missing.');
    end
end
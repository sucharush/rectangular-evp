function [lam_star, sig_star, info] = minimizer_fminsearch(a, b, objective, opts)
% Local minimization on [a,b] using fminsearch, initialized at midpoint.
%
% Required:
%   a, b          : scalar bracket endpoints
%   objective.eval: function handle
%   objective.get_count: function handle
%
% Optional in opts:
%   tol_x
%   tol_fun
%   max_iter
%   max_fun_evals

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts, 'tol_x'),         opts.tol_x = 1e-13; end
    if ~isfield(opts, 'tol_fun'),       opts.tol_fun = 1e-12; end
    if ~isfield(opts, 'max_iter'),      opts.max_iter = 200; end
    if ~isfield(opts, 'max_fun_evals'), opts.max_fun_evals = 500; end

    x0 = 0.5 * (a + b);

    options = optimset( ...
        'TolX', opts.tol_x, ...
        'TolFun', opts.tol_fun, ...
        'MaxIter', opts.max_iter, ...
        'MaxFunEvals', opts.max_fun_evals, ...
        'Display', 'off');

    lam_star = NaN;
    sig_star = Inf;
    status = 'ok';
    message = '';

    try
        [lam_star, sig_star] = fminsearch(objective.eval, x0, options);
    catch ME
        status = 'exception';
        message = ME.message;
    end

    info = struct();
    info.method = 'fminsearch';
    info.status = status;
    info.message = message;
    info.bracket = [a, b];
    info.x0 = x0;
    info.n_calls = objective.get_count();
end
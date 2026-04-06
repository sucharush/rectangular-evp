function [lam_star, sig_star, info] = minimizer_fminbnd(a, b, objective, opts)
% Local bounded minimization on [a,b] using fminbnd.

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts, 'tol_x'),    opts.tol_x = 1e-13; end
    if ~isfield(opts, 'max_iter'), opts.max_iter = 200; end

    options = optimset( ...
        'TolX', opts.tol_x, ...
        'MaxIter', opts.max_iter, ...
        'Display', 'off');

    lam_star = NaN;
    sig_star = Inf;
    status = 'ok';
    message = '';

    try
        [lam_star, sig_star] = fminbnd(objective.eval, a, b, options);
    catch ME
        status = 'exception';
        message = ME.message;
    end

    info = struct();
    info.method = 'fminbnd';
    info.status = status;
    info.message = message;
    info.bracket = [a, b];
    info.x0 = [];
    info.n_calls = objective.get_count();
end
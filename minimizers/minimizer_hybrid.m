function [lam_star, sig_star, info] = minimizer_hybrid(a, b, objective, opts)
% Hybrid local minimizer:
% 1. a few trisection steps to shrink bracket
% 2. fminsearch from midpoint of reduced interval

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts, 'tol_x'),         opts.tol_x = 1e-13; end
    if ~isfield(opts, 'tol_fun'),       opts.tol_fun = 1e-12; end
    if ~isfield(opts, 'max_iter'),      opts.max_iter = 200; end
    if ~isfield(opts, 'max_fun_evals'), opts.max_fun_evals = 500; end
    if ~isfield(opts, 'n_pre'),         opts.n_pre = 4; end

    a0 = a;
    b0 = b;

    status = 'ok';
    message = '';

    try
        for k = 1:opts.n_pre
            if b - a <= opts.tol_x
                break;
            end

            x1 = a + (b - a) / 3;
            x2 = b - (b - a) / 3;

            f1 = objective.eval(x1);
            f2 = objective.eval(x2);

            if f1 <= f2
                b = x2;
            else
                a = x1;
            end
        end

        x0 = 0.5 * (a + b);

        options = optimset( ...
            'TolX', opts.tol_x, ...
            'TolFun', opts.tol_fun, ...
            'MaxIter', opts.max_iter, ...
            'MaxFunEvals', opts.max_fun_evals, ...
            'Display', 'off');

        [lam_star, sig_star] = fminsearch(objective.eval, x0, options);

    catch ME
        lam_star = NaN;
        sig_star = Inf;
        x0 = NaN;
        status = 'exception';
        message = ME.message;
    end

    info = struct();
    info.method = 'hybrid';
    info.status = status;
    info.message = message;
    info.original_bracket = [a0, b0];
    info.reduced_bracket = [a, b];
    info.x0 = x0;
    info.n_calls = objective.get_count();
end
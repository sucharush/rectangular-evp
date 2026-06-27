function [lam_star, sig_star, info] = minimizer_golden_section(a, b, objective, opts)
% Golden-section search for a unimodal objective on [a,b].

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts, 'tol_x'),    opts.tol_x = 1e-13; end
    if ~isfield(opts, 'max_iter'), opts.max_iter = 200; end

    status = 'ok';
    message = '';

    try
        r = (sqrt(5) - 1) / 2;

        x1 = b - r * (b - a);
        x2 = a + r * (b - a);
        f1 = objective.eval(x1);
        f2 = objective.eval(x2);

        for k = 1:opts.max_iter
            if b - a <= opts.tol_x
                break;
            end

            if f1 <= f2
                b = x2;
                x2 = x1;
                f2 = f1;

                x1 = b - r * (b - a);
                f1 = objective.eval(x1);
            else
                a = x1;

                x1 = x2;
                f1 = f2;

                x2 = a + r * (b - a);
                f2 = objective.eval(x2);
            end
        end

        lam_star = 0.5 * (a + b);
        sig_star = objective.eval(lam_star);

    catch ME
        lam_star = NaN;
        sig_star = Inf;
        status = 'exception';
        message = ME.message;
    end

    info = struct();
    info.method = 'golden_section';
    info.status = status;
    info.message = message;
    info.bracket = [a, b];
    info.x0 = [];
    info.n_calls = objective.get_count();
end

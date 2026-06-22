function [lam_star, sig_star, info] = minimizer_aaa_ellipse(a, b, objective, opts)
% AAA-based local minimizer using random Bernstein ellipse samples.

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts, 'nZ'),       opts.nZ = 100; end
    if ~isfield(opts, 'rho'),      opts.rho = 1.05; end
    if ~isfield(opts, 'delta'),    opts.delta = 1e-13; end
    if ~isfield(opts, 'imag_tol'), opts.imag_tol = 1e-3; end
    if ~isfield(opts, 'mmax'),     opts.mmax = 100; end
    if ~isfield(opts, 'eval_tol'), opts.eval_tol = 1e-13; end

    lam_star = NaN;
    sig_star = Inf;
    status = 'ok';
    message = '';
    candidates = [];

    try
        if exist('aaa', 'file') ~= 2
            error('minimizer_aaa_ellipse:MissingAAA', ...
                'aaa.m is not on the MATLAB path.');
        end

        Z = local_bernstein_ellipse_points(a, b, opts.rho, opts.nZ);
        F = arrayfun(@(z) objective.eval(z), Z);
        Rdata = 1 ./ (F + opts.delta);

        [~, pol] = aaa(Rdata, Z, 'mmax', opts.mmax, 'tol', opts.eval_tol);
        candidates = local_real_candidates(pol, a, b, opts.imag_tol);

        if isempty(candidates)
            status = 'no_candidate';
            message = 'AAA returned no near-real pole in the bracket.';
        else
            s_at = arrayfun(@(x) objective.eval(x), candidates);
            [sig_star, idx_best] = min(s_at);
            lam_star = candidates(idx_best);
        end
    catch ME
        status = 'exception';
        message = ME.message;
    end

    info = struct();
    info.method = 'aaa_ellipse';
    info.status = status;
    info.message = message;
    info.bracket = [a, b];
    info.x0 = [];
    info.n_calls = objective.get_count();
    info.candidates = candidates;
end


% function Z = local_bernstein_ellipse_points(a, b, rho, n)
%     theta = 2*pi .* rand(n+1, 1).';
%     theta(end) = [];
% 
%     w = rho * exp(1i * theta);
%     xi = 0.5 * (w + 1 ./ w);
%     Z = 0.5 * (a + b) + 0.5 * (b - a) * xi;
% end

function Z = local_bernstein_ellipse_points(a, b, rho, n)
    A = 0.5 * (rho + 1/rho);
    B = 0.5 * (rho - 1/rho);

    r_rand = sqrt(rand(1, n)); 
    theta = 2*pi * rand(1, n);
    
    X = A * r_rand .* cos(theta);
    Y = B * r_rand .* sin(theta);
    xi = X + 1i * Y;
    
    Z = 0.5*(a+b) + 0.5*(b-a)*xi;
end


function candidates = local_real_candidates(pol, a, b, imag_tol)
    candidates = pol(abs(imag(pol)) < imag_tol);
    candidates = real(candidates);
    candidates = candidates(candidates >= a & candidates <= b);
    candidates = sort(candidates(:).');
end

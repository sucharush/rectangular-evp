function scan = scan_sigma(sigma_fun, opts)
% Global scan of sigma(lambda) on a user-supplied grid.
%
% Required:
%   sigma_fun (function handle)
%   opts.lamvec
%
% Optional:
%   opts.detect_mode = 'strict_local_min'   (default)
%
% Output:
%   scan.lamvec
%   scan.S
%   scan.candidate_idx

    if ~isfield(opts, 'lamvec')
        error('scan_sigma: opts.lamvec is required.');
    end
    if ~isfield(opts, 'detect_mode')
        opts.detect_mode = 'strict_local_min';
    end

    lamvec = opts.lamvec(:).';   % keep row vector style
    S = zeros(size(lamvec));

    for i = 1:numel(lamvec)
        % S(i) = problem.ops.sigma(lamvec(i));
        S(i) = sigma_fun(lamvec(i));
    end

    switch lower(opts.detect_mode)
        case 'strict_local_min'
            J = 2:numel(lamvec)-1;
            J = J(S(J) < S(J-1) & S(J) < S(J+1));

        otherwise
            error('scan_sigma: unknown detect_mode "%s".', opts.detect_mode);
    end

    scan = struct();
    scan.lamvec = lamvec;
    scan.S = S;
    scan.candidate_idx = J;
end
% function out = aaa_rect_pencil_interval(F_op, a, b, opts)
% % Build a Theorem-3-style linear surrogate pencil for a matrix-valued
% % operator F(lambda) on a bounded interval [a,b].
% %
% % INPUT
% %   F_op : function handle, F_op(lam) returns a p-by-q matrix
% %   a,b  : interval endpoints
% %   opts : struct with AAA / sketch settings
% %
% % OUTPUT
% %   out : struct containing support points, weights, Newton parameters,
% %         sampled matrices, and the surrogate pencil A - lam B
% 
%     if nargin < 4
%         opts = struct();
%     end
% 
%     % defaults
%     if ~isfield(opts, 'nCand'),      opts.nCand = 400; end
%     if ~isfield(opts, 'aaa_tol'),    opts.aaa_tol = 1e-12; end
%     if ~isfield(opts, 'mmax'),       opts.mmax = 50; end
%     if ~isfield(opts, 'seed'),       opts.seed = 1; end
%     if ~isfield(opts, 'ell'),        opts.ell = 2; end
%     if ~isfield(opts, 'max_norm'),   opts.max_norm = true; end
%     if ~isfield(opts, 'svd_update'), opts.svd_update = true; end
%     if ~isfield(opts, 'beta_mode'),  opts.beta_mode = 'ones'; end
%     if ~isfield(opts, 'rho'),        opts.rho = 1.05; end
%     if ~isfield(opts, 'sample_mode'), opts.sample_mode = 'real'; end
% 
%     % 0. size check
%     Fmid = F_op((a+b)/2);
%     [p, q] = size(Fmid);
% 
%     % 1. random sketch probes
%     rng(opts.seed);
% 
%     ell = opts.ell;
%     U = randn(p, ell) + 1i * randn(p, ell);
%     V = randn(q, ell) + 1i * randn(q, ell);
% 
%     for i = 1:ell
%         U(:,i) = U(:,i) / norm(U(:,i));
%         V(:,i) = V(:,i) / norm(V(:,i));
%     end
% 
%     % 2. candidate sample points
%     switch lower(opts.sample_mode)
%         case 'real'
%             Z = linspace(a, b, opts.nCand).';
% 
%         case 'ellipse'
%             Z = aaa_bernstein_ellipse_points(a, b, opts.rho, opts.nCand).';
% 
%         otherwise
%             error('aaa_rect_pencil_interval: unknown sample_mode "%s".', opts.sample_mode);
%     end
% 
%     % 3. build sketched scalar surrogates
%     GZ = aaa_build_sketch_surrogates(F_op, Z, U, V);
% 
%     scale_factor = max(abs(GZ(:)));
%     if scale_factor == 0
%         scale_factor = 1;
%     end
%     GZ_normalized = GZ / scale_factor;
% 
%     aaa_opts = struct();
%     aaa_opts.tol = opts.aaa_tol;
%     aaa_opts.max_norm = opts.max_norm;
%     aaa_opts.svd_update = opts.svd_update;
% 
%     % IMPORTANT:
%     % sketchAAA must be on path
%     [zj, wj, ind, stats_aaa, statsT_aaa] = sketchAAA(GZ_normalized, Z, opts.mmax, aaa_opts);
% 
%     zj  = zj(:);
%     wj  = wj(:);
%     ind = ind(:);
% 
%     m = numel(zj) - 1;
%     if m < 1
%         error('aaa_rect_pencil_interval: sketchAAA returned degree m < 1.');
%     end
% 
%     % 4. barycentric -> Newton parameters
%     sigma = zj(1:m);
%     beta  = zeros(m,1);
%     h     = zeros(m,1);
%     k     = zeros(m,1);
%     xi    = zeros(m,1);
% 
%     for j = 1:m
%         ratio = -wj(j) / wj(j+1);
% 
%         switch lower(opts.beta_mode)
%             case 'ones'
%                 beta(j) = 1;
%                 k(j) = ratio;
%                 h(j) = zj(j+1) * ratio;
% 
%             case 'unitk'
%                 k(j) = 1;
%                 beta(j) = ratio;
%                 h(j) = zj(j+1);
% 
%             otherwise
%                 error('aaa_rect_pencil_interval: unknown beta_mode "%s".', opts.beta_mode);
%         end
% 
%         xi(j) = h(j) / k(j);
%     end
% 
%     % 5. sampled operator blocks
%     D = cell(m+1,1);
%     for j = 1:m+1
%         Dj = F_op(zj(j));
%         [pj, qj] = size(Dj);
% 
%         if pj ~= p || qj ~= q
%             error('aaa_rect_pencil_interval: inconsistent operator size across samples.');
%         end
% 
%         D{j} = Dj;
%     end
% 
%     % 6. build surrogate pencil
%     [A, B] = build_rectangular_linear_pencil(D, sigma, beta, h, k);
% 
%     % 7. pack output
%     out = struct();
%     out.a = a;
%     out.b = b;
% 
%     out.p = p;
%     out.q = q;
%     out.m = m;
% 
%     out.ell = ell;
%     out.U = U;
%     out.V = V;
% 
%     out.Z = Z;
%     out.GZ = GZ;
%     out.GZ_normalized = GZ_normalized;
%     out.scale_factor = scale_factor;
%     out.ind = ind;
% 
%     out.zj = zj;
%     out.wj = wj;
% 
%     out.sigma = sigma;
%     out.beta = beta;
%     out.h = h;
%     out.k = k;
%     out.xi = xi;
% 
%     out.D = D;
%     out.A = A;
%     out.B = B;
% 
%     out.stats_aaa = stats_aaa;
%     out.statsT_aaa = statsT_aaa;
% end
% 


function out = aaa_rect_pencil_interval(op, a, b, opts)
% Unified AAA builder for rectangular matrix-valued operators.
%
% INPUT
%   op   : function handle
%   a,b  : interval endpoints
%   opts : struct
%
% Required mode:
%   opts.method = 'direct'
%       op is a matrix-valued operator QB_op(lambda)
%
%   opts.method = 'procrustes_from_A'
%       op is A_op(lambda); routine builds QR factors, fixes signs,
%       aligns them by Procrustes, and uses the aligned top block.
%       Requires opts.mB.
%
% Common optional fields:
%   opts.nCand
%   opts.aaa_tol
%   opts.mmax
%   opts.seed
%   opts.ell
%   opts.max_norm
%   opts.svd_update
%   opts.beta_mode
%   opts.sample_mode = 'real' or 'ellipse'
%   opts.rho
%
% Additional optional fields for 'procrustes_from_A':
%   opts.sign_fix = true/false
%   opts.store_Q = true/false
%
% OUTPUT
%   out.method
%   out.Z
%   out.D
%   out.zj, out.wj
%   out.sigma, out.beta, out.h, out.k, out.xi
%   out.A, out.B
%   and, for procrustes mode:
%   out.Q_all
%   out.QB_all

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts, 'method')
        opts.method = 'direct';
    end
    if ~isfield(opts, 'nCand'),      opts.nCand = 400; end
    if ~isfield(opts, 'aaa_tol'),    opts.aaa_tol = 1e-12; end
    if ~isfield(opts, 'mmax'),       opts.mmax = 50; end
    if ~isfield(opts, 'seed'),       opts.seed = 1; end
    if ~isfield(opts, 'ell'),        opts.ell = 2; end
    if ~isfield(opts, 'max_norm'),   opts.max_norm = true; end
    if ~isfield(opts, 'svd_update'), opts.svd_update = true; end
    if ~isfield(opts, 'beta_mode'),  opts.beta_mode = 'ones'; end
    if ~isfield(opts, 'sample_mode'), opts.sample_mode = 'real'; end
    if ~isfield(opts, 'rho'),        opts.rho = 1.05; end
    if ~isfield(opts, 'sign_fix'),   opts.sign_fix = true; end
    if ~isfield(opts, 'store_Q'),    opts.store_Q = true; end

    % ------------------------------------------------------------
    % 0. candidate sample points
    % ------------------------------------------------------------
    switch lower(opts.sample_mode)
        case 'real'
            Z = linspace(a, b, opts.nCand).';
        case 'ellipse'
            Z = aaa_bernstein_ellipse_points(a, b, opts.rho, opts.nCand).';
        otherwise
            error('aaa_rect_pencil_interval: unknown sample_mode "%s".', opts.sample_mode);
    end

    % ------------------------------------------------------------
    % 1. build sampled operator blocks D_j
    % ------------------------------------------------------------
    switch lower(opts.method)
        case 'direct'
            [D_all, aux] = local_sample_direct(op, Z);

        case 'procrustes_from_a'
            if ~isfield(opts, 'mB') || isempty(opts.mB)
                error('aaa_rect_pencil_interval: opts.mB is required for method="procrustes_from_A".');
            end
            [D_all, aux] = local_sample_procrustes_from_A(op, Z, opts);

        otherwise
            error('aaa_rect_pencil_interval: unknown method "%s".', opts.method);
    end

    [p, q, nZ] = size(D_all);

    if p < q
        warning('aaa_rect_pencil_interval: sampled operator appears underdetermined (p < q).');
    end

    % ------------------------------------------------------------
    % 2. build sketch surrogates
    % ------------------------------------------------------------
    rng(opts.seed);

    ell = opts.ell;
    U = randn(p, ell) + 1i * randn(p, ell);
    V = randn(q, ell) + 1i * randn(q, ell);

    for i = 1:ell
        U(:,i) = U(:,i) / norm(U(:,i));
        V(:,i) = V(:,i) / norm(V(:,i));
    end

    GZ = zeros(ell, nZ);
    for j = 1:nZ
        Fj = D_all(:,:,j);
        for i = 1:ell
            GZ(i,j) = U(:,i)' * Fj * V(:,i);
        end
    end

    scale_factor = max(abs(GZ(:)));
    if scale_factor == 0
        scale_factor = 1;
    end
    GZ_normalized = GZ / scale_factor;

    aaa_opts = struct();
    aaa_opts.tol = opts.aaa_tol;
    aaa_opts.max_norm = opts.max_norm;
    aaa_opts.svd_update = opts.svd_update;

    [zj, wj, ind, stats_aaa, statsT_aaa] = sketchAAA(GZ_normalized, Z, opts.mmax, aaa_opts);

    zj  = zj(:);
    wj  = wj(:);
    ind = ind(:);

    m = numel(zj) - 1;
    if m < 1
        error('aaa_rect_pencil_interval: sketchAAA returned degree m < 1.');
    end

    % ------------------------------------------------------------
    % 3. barycentric -> Newton parameters
    % ------------------------------------------------------------
    sigma = zj(1:m);
    beta  = zeros(m,1);
    h     = zeros(m,1);
    k     = zeros(m,1);
    xi    = zeros(m,1);

    for j = 1:m
        ratio = -wj(j) / wj(j+1);

        switch lower(opts.beta_mode)
            case 'ones'
                beta(j) = 1;
                k(j) = ratio;
                h(j) = zj(j+1) * ratio;

            case 'unitk'
                k(j) = 1;
                beta(j) = ratio;
                h(j) = zj(j+1);

            otherwise
                error('aaa_rect_pencil_interval: unknown beta_mode "%s".', opts.beta_mode);
        end

        xi(j) = h(j) / k(j);
    end

    % ------------------------------------------------------------
    % 4. support-point matrix blocks
    % ------------------------------------------------------------
    D = cell(m+1,1);
    for j = 1:m+1
        D{j} = D_all(:,:,ind(j));
    end

    % ------------------------------------------------------------
    % 5. surrogate pencil
    % ------------------------------------------------------------
    [A, B] = build_rectangular_linear_pencil(D, sigma, beta, h, k);

    % ------------------------------------------------------------
    % 6. output
    % ------------------------------------------------------------
    out = struct();
    out.method = lower(opts.method);

    out.a = a;
    out.b = b;
    out.p = p;
    out.q = q;
    out.m = m;

    out.ell = ell;
    out.U = U;
    out.V = V;

    out.Z = Z;
    out.GZ = GZ;
    out.GZ_normalized = GZ_normalized;
    out.ind = ind;
    out.scale_factor = scale_factor;

    out.zj = zj;
    out.wj = wj;

    out.sigma = sigma;
    out.beta = beta;
    out.h = h;
    out.k = k;
    out.xi = xi;

    out.D = D;
    out.A = A;
    out.B = B;

    out.stats_aaa = stats_aaa;
    out.statsT_aaa = statsT_aaa;

    % mode-specific stored data
    if strcmpi(out.method, 'procrustes_from_a')
        out.mB = opts.mB;
        out.sign_fix = opts.sign_fix;
        out.QB_all = aux.QB_all;
        if isfield(aux, 'Q_all')
            out.Q_all = aux.Q_all;
        end
    end
end


function [D_all, aux] = local_sample_direct(QB_op, Z)
    F1 = QB_op(Z(1));
    [p, q] = size(F1);
    nZ = numel(Z);

    D_all = zeros(p, q, nZ, class(F1));
    D_all(:,:,1) = F1;

    for j = 2:nZ
        Fj = QB_op(Z(j));
        [pj, qj] = size(Fj);
        if pj ~= p || qj ~= q
            error('local_sample_direct: inconsistent matrix size across samples.');
        end
        D_all(:,:,j) = Fj;
    end

    aux = struct();
end


function [D_all, aux] = local_sample_procrustes_from_A(A_op, Z, opts)
    A1 = A_op(Z(1));
    [Q1, R1] = qr(A1, 0);

    if opts.sign_fix
        d = sign(diag(R1));
        d(d == 0) = 1;
        Q1 = Q1 * diag(d);
    end

    QB1 = Q1(1:opts.mB, :);
    [p, q] = size(QB1);
    nZ = numel(Z);

    D_all = zeros(p, q, nZ, class(QB1));
    D_all(:,:,1) = QB1;

    if opts.store_Q
        Q_all = zeros(size(Q1,1), size(Q1,2), nZ, class(Q1));
        Q_all(:,:,1) = Q1;
    end

    Q_prev = Q1;

    for j = 2:nZ
        Aj = A_op(Z(j));
        [Qj, Rj] = qr(Aj, 0);

        if opts.sign_fix
            d = sign(diag(Rj));
            d(d == 0) = 1;
            Qj = Qj * diag(d);
        end

        [U, ~, V] = svd(Q_prev' * Qj, 'econ');
        Ralign = V * U';
        Qj = Qj * Ralign;

        D_all(:,:,j) = Qj(1:opts.mB, :);

        if opts.store_Q
            Q_all(:,:,j) = Qj;
        end

        Q_prev = Qj;
    end

    aux = struct();
    aux.QB_all = D_all;
    if opts.store_Q
        aux.Q_all = Q_all;
    end
end

function [A, B] = build_rectangular_linear_pencil(D, sigma, beta, h, k)
% Direct rectangular analogue of Theorem 3.
%
% D{1},...,D{m+1} correspond to D_0,...,D_m, each p-by-q
% sigma, beta, h, k are length-m

    m = numel(beta);
    [p,q] = size(D{1});

    nrows = p + (m-1)*q;
    ncols = m*q;

    A = sparse(nrows, ncols);
    B = sparse(nrows, ncols);

    hm = h(m);
    km = k(m);
    betam = beta(m);

    % top block row
    for j = 1:(m-1)
        cols = (j-1)*q + (1:q);
        A(1:p, cols) = hm * D{j};
        B(1:p, cols) = km * D{j};
    end

    cols = (m-1)*q + (1:q);
    % FIX 1: Removed spurious 'hm' from the D{m+1} tail correction
    A(1:p, cols) = hm * D{m} - (sigma(m) / betam) * D{m+1};
    B(1:p, cols) = km * D{m} - (1 / betam) * D{m+1};

    Iq = speye(q);

    % lower block rows
    for i = 1:(m-1)
        rows      = p + (i-1)*q + (1:q);
        col_left  = (i-1)*q + (1:q);
        col_right = i*q     + (1:q);

        % no 'h(i)' from the last term
        A(rows, col_left)  = sigma(i) * Iq;
        A(rows, col_right) = h(i) * beta(i)  * Iq;

        B(rows, col_left)  = 1 * Iq;
        B(rows, col_right) = k(i) * beta(i)  * Iq;
    end
    % Compute a rough norm of the top row vs the identity blocks
    norm_top = norm(A(1:p, :), 'inf') + norm(B(1:p, :), 'inf');
    norm_bot = norm(A(p+1:end, :), 'inf') + norm(B(p+1:end, :), 'inf');

    % Scale the top row equations to match the lower recurrence equations
    gamma = norm_bot / max(norm_top, 1e-14);
    A(1:p, :) = gamma * A(1:p, :);
    B(1:p, :) = gamma * B(1:p, :);
end
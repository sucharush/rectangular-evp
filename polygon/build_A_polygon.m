function A = build_A_polygon(P, corner_list, Mcorner, lam)
% Build the oversampled collocation matrix A(lambda)

    Z = sqrt(lam);
    npts = size(P,1);
    blocks = cell(numel(corner_list),1);

    for j = 1:numel(corner_list)
        sj = corner_list{j};
        c = sj.c;
        e0 = sj.e0;
        alpha = sj.alpha;

        d = P - c;
        rho = sqrt(sum(d.^2, 2));

        detv = e0(1)*d(:,2) - e0(2)*d(:,1);
        dotv = e0(1)*d(:,1) + e0(2)*d(:,2);

        phi0 = atan2(detv, dotv);
        cands = [phi0 - 2*pi, phi0, phi0 + 2*pi];

        dist = max(0, -cands) + max(0, cands - sj.omega);
        [~, idx] = min(dist, [], 2);
        phi = cands(sub2ind(size(cands), (1:npts)', idx));

        m = 1:Mcorner;
        nu = alpha * m;

        Rmat  = (Z * rho) * ones(1, Mcorner);
        NUmat = ones(npts,1) * nu;

        blocks{j} = besselj(NUmat, Rmat) .* sin(NUmat .* (phi * ones(1, Mcorner)));
    end

    A = [blocks{:}];
end
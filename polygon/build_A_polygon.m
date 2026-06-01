function A = build_A_polygon(geom, lam)
% Build the oversampled collocation matrix A(lambda)

    Z = sqrt(lam);
    Mcorner = geom.Mcorner;
    corner_geom = geom.corner_geom;
    npts = size(geom.P, 1);
    blocks = cell(numel(corner_geom), 1);
    m = 1:Mcorner;

    for j = 1:numel(corner_geom)
        gj = corner_geom{j};
        nu = gj.alpha * m;

        Rmat  = (Z * gj.rho) * ones(1, Mcorner);
        NUmat = ones(npts,1) * nu;

        blocks{j} = besselj(NUmat, Rmat) .* sin(NUmat .* (gj.phi * ones(1, Mcorner)));
    end

    A = [blocks{:}];
end

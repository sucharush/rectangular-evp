function A = build_A_polygon(geom, lam)
% Build the oversampled collocation matrix A(lambda)

    Z = sqrt(lam);
    corner_geom = geom.corner_geom;
    blocks = cell(numel(corner_geom), 1);

    for j = 1:numel(corner_geom)
        gj = corner_geom{j};
        Rmat = (Z * gj.rho) * ones(1, size(gj.nu_matrix, 2));
        blocks{j} = besselj(gj.nu_matrix, Rmat) .* gj.sin_term;
    end

    A = [blocks{:}];
end

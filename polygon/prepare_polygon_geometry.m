function geom = prepare_polygon_geometry(V, cfg)
    if any(V(1,:) ~= V(end,:))
        V = [V; V(1,:)];
    end

    area2 = sum(V(1:end-1,1).*V(2:end,2) - V(2:end,1).*V(1:end-1,2));
    if area2 < 0
        V = flipud(V);
    end

    B = cfg.sampling.boundary_fun(V, cfg.nb_per_edge, cfg);
    I = cfg.sampling.interior_fun(V, cfg.nI, cfg);
    P = [B; I];

    corner_list = detect_singular_corners(V);
    corner_geom = prepare_corner_geometry(P, corner_list, cfg.Mcorner);

    geom = struct();
    geom.V = V;
    geom.B = B;
    geom.I = I;
    geom.mB = size(B,1);
    geom.corner_list = corner_list;
    geom.corner_geom = corner_geom;
end


function corner_geom = prepare_corner_geometry(P, corner_list, Mcorner)
    npts = size(P, 1);
    corner_geom = cell(numel(corner_list), 1);
    m = 1:Mcorner;

    for j = 1:numel(corner_list)
        sj = corner_list{j};
        d = P - sj.c;
        rho = sqrt(sum(d.^2, 2));

        detv = sj.e0(1) * d(:,2) - sj.e0(2) * d(:,1);
        dotv = sj.e0(1) * d(:,1) + sj.e0(2) * d(:,2);

        phi0 = atan2(detv, dotv);
        cands = [phi0 - 2*pi, phi0, phi0 + 2*pi];

        dist = max(0, -cands) + max(0, cands - sj.omega);
        [~, idx] = min(dist, [], 2);
        phi = cands(sub2ind(size(cands), (1:npts)', idx));

        nu_matrix = ones(npts, 1) * (sj.alpha * m);
        sin_term = sin(nu_matrix .* (phi * ones(1, Mcorner)));

        corner_geom{j} = struct( ...
            'rho', rho, ...
            'nu_matrix', nu_matrix, ...
            'sin_term', sin_term);
    end
end

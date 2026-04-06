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

    geom = struct();
    geom.V = V;
    geom.B = B;
    geom.I = I;
    geom.P = P;
    geom.mB = size(B,1);
    geom.corner_list = corner_list;
    geom.Mcorner = cfg.Mcorner;
end
function problem = build_polygon_problem(cfg)
    geom = prepare_polygon_geometry(cfg.V, cfg);

    problem = struct();
    problem.name = cfg.name;
    problem.family = 'polygon_laplace';

    problem.data = struct();
    problem.data.cfg = cfg;
    problem.data.geom = geom;

    problem.ops = struct();
    problem.ops.A = @(lam) build_A_polygon(geom.P, geom.corner_list, cfg.Mcorner, lam);
    problem.ops.QB = @(lam) build_QB_polygon(geom, cfg, lam);
    problem.ops.sigma = @(lam) sigma_polygon(geom, cfg, lam);
    % general solver-facing matrix-valued operator
    problem.ops.F = problem.ops.QB;
    
    problem.meta = struct();
    problem.meta.mB = geom.mB;
    problem.meta.nI = size(geom.I, 1);
    problem.meta.nCorners = numel(geom.corner_list);
    problem.meta.nominal_basis_dim = cfg.Mcorner * numel(geom.corner_list);

end
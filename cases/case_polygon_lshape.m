function cfg = case_polygon_lshape()
    cfg = struct();
    cfg.name = 'polygon_lshape';

    % closed polygon
    cfg.V = [ ...
         0  0;
         0  1;
        -1  1;
        -1 -1;
         1 -1;
         1  0;
         0  0];

    cfg.Mcorner = 40;
    cfg.nb_per_edge = 40;
    cfg.nI = 50;

    % sampling policy
    cfg.sampling = struct();
    cfg.sampling.boundary_fun = @sample_boundary_chebyshev;
    cfg.sampling.interior_fun = @sample_interior_uniform;

    % operator / QR policy
    cfg.normalize_columns = true;
    cfg.qr_tau = 1e-13;
end
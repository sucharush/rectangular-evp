function cfg = case_polygon_hshape()
    a = 2;
    l = 1;
    t = 0.3;

    cfg = struct();
    cfg.name = 'polygon_hshape';

    cfg.V = [ ...
      0,0;
      a,0;
      a,a/2-t/2;
      a+l,a/2-t/2;
      a+l,0;
      2*a+l,0;
      2*a+l,a;
      a+l,a;
      a+l,a/2+t/2;
      a,a/2+t/2;
      a,a;
      0,a;
      0,0];

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
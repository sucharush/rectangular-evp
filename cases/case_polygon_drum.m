function cfg = case_polygon_drum(side)
    if nargin < 1
        side = 'left';
    end
    cfg = struct();
    cfg.name = 'polygon_drum';

    V1 = [ ...
        -1  1
        -1  3
        -3  1
        -1 -1
         1 -1
         1 -3
         3 -1
         3  1
        -1  1 ];
    
    V2 = [ ...
        -1  1
        -1  3
        -3  3
        -3  1
         1 -3
         1 -1
         3 -1
         1  1
        -1  1 ];
    switch lower(side)
        case 'left'
            cfg.V = V1;
        case 'right'
            cfg.V = V2;
        otherwise
            cfg.V = V1;
            print('undefined type, use ''left'' by default');
    end
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
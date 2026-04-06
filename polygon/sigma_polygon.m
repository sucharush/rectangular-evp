function s = sigma_polygon(geom, cfg, lam)
    QB = build_QB_polygon(geom, cfg, lam);

    if isempty(QB)
        s = NaN;
        return;
    end

    s = min(svd(QB, 'econ'));
end
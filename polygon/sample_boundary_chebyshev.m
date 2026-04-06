function B = sample_boundary_chebyshev(V, nb_per_edge, cfg) %#ok<INUSD>
% V is assumed closed: V(end,:) = V(1,:)

    n = nb_per_edge;
    j = (1:n)';
    s = 0.5 * (1 - cos((2*j - 1) * pi / (2*n)));  % in (0,1)

    B = [];

    for e = 1:size(V,1)-1
        a = V(e,:);
        b = V(e+1,:);
        pts = a + s .* (b - a);
        B = [B; pts]; %#ok<AGROW>
    end
end
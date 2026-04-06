function Rz = aaa_eval_matrix_barycentric(z, zj, wj, D)
    hit = find(abs(z - zj) < 1e-14 * max(1, abs(z)), 1);
    if ~isempty(hit)
        Rz = D{hit};
        return;
    end

    [p, q] = size(D{1});
    num = zeros(p, q);
    den = 0;

    for j = 1:numel(zj)
        alpha = wj(j) / (z - zj(j));
        num = num + alpha * D{j};
        den = den + alpha;
    end

    Rz = num / den;
end
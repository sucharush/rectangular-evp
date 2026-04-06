function Rz = aaa_eval_matrix_newton(z, D, sigma, beta, h, k)
% D{1},...,D{m+1} correspond to D_0,...,D_m

    m = numel(beta);

    if numel(sigma) ~= m || numel(h) ~= m || numel(k) ~= m
        error('aaa_eval_matrix_newton: sigma, beta, h, k must all have length m.');
    end
    if numel(D) ~= m + 1
        error('aaa_eval_matrix_newton: D must contain m+1 blocks.');
    end

    b = zeros(m+1, 1);
    b(1) = 1;

    for j = 1:m
        b(j+1) = ((z - sigma(j)) / (beta(j) * (h(j) - k(j) * z))) * b(j);
    end

    [p, q] = size(D{1});
    Rz = zeros(p, q);

    for j = 1:m+1
        Rz = Rz + b(j) * D{j};
    end
end
function Z = aaa_bernstein_ellipse_points(a, b, rho, n)
    theta = 2*pi .* rand(n+1, 1).';
    theta(end) = [];

    w = rho * exp(1i * theta);
    xi = 0.5 * (w + 1 ./ w);
    Z = 0.5 * (a + b) + 0.5 * (b - a) * xi;
end
function I = sample_interior_uniform(V, nI, cfg) %#ok<INUSD>
% Rejection sampling in bounding box

    xmin = min(V(:,1));
    xmax = max(V(:,1));
    ymin = min(V(:,2));
    ymax = max(V(:,2));

    I = zeros(nI, 2);
    cnt = 0;

    while cnt < nI
        cand = [ ...
            xmin + (xmax - xmin) * rand(5000,1), ...
            ymin + (ymax - ymin) * rand(5000,1)];

        in = inpolygon(cand(:,1), cand(:,2), V(:,1), V(:,2));
        cand = cand(in,:);

        take = min(size(cand,1), nI - cnt);
        I(cnt+1:cnt+take, :) = cand(1:take, :);
        cnt = cnt + take;
    end
end
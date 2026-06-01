function I = sample_interior_uniform(V, nI, cfg)
% Rejection sampling in bounding box, with optional targeted box sampling.

    xmin = min(V(:,1));
    xmax = max(V(:,1));
    ymin = min(V(:,2));
    ymax = max(V(:,2));

    n_target = local_target_count(cfg, nI);
    n_global = nI - n_target;

    I = zeros(nI, 2);
    cnt = 0;

    while cnt < n_global
        cand = [ ...
            xmin + (xmax - xmin) * rand(5000,1), ...
            ymin + (ymax - ymin) * rand(5000,1)];

        in = inpolygon(cand(:,1), cand(:,2), V(:,1), V(:,2));
        cand = cand(in,:);

        take = min(size(cand,1), n_global - cnt);
        I(cnt+1:cnt+take, :) = cand(1:take, :);
        cnt = cnt + take;
    end

    if n_target > 0
        box = cfg.sampling.channel_box;
        while cnt < nI
            cand = [ ...
                box(1) + (box(2) - box(1)) * rand(5000,1), ...
                box(3) + (box(4) - box(3)) * rand(5000,1)];

            in = inpolygon(cand(:,1), cand(:,2), V(:,1), V(:,2));
            cand = cand(in,:);

            take = min(size(cand,1), nI - cnt);
            I(cnt+1:cnt+take, :) = cand(1:take, :);
            cnt = cnt + take;
        end
    end
end


function n_target = local_target_count(cfg, nI)
    n_target = 0;
    if ~isfield(cfg, 'sampling') || ~isfield(cfg.sampling, 'channel_box')
        return;
    end

    frac = 0.1;
    if isfield(cfg.sampling, 'frac_channel') && ~isempty(cfg.sampling.frac_channel)
        frac = cfg.sampling.frac_channel;
    end

    frac = max(0, min(1, frac));
    n_target = round(frac * nI);
end

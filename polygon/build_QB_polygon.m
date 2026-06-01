function QB = build_QB_polygon(geom, cfg, lam)
% Centralized place for:
% 1. build A(lambda)
% 2. optional column normalization
% 3. QR truncation
% 4. return boundary block QB

    A = build_A_polygon(geom, lam);
    qr_opts = struct( ...
        'normalize_columns', is_true_field(cfg, 'normalize_columns'), ...
        'pivot', is_true_field(cfg, 'qr_pivot'), ...
        'sign_fix', is_true_field(cfg, 'qr_sign_fix'), ...
        'qr_tau', get_field_or_empty(cfg, 'qr_tau'));

    QB = build_QB_from_A(A, geom.mB, qr_opts);
end


function tf = is_true_field(s, name)
    tf = isfield(s, name) && ~isempty(s.(name)) && s.(name);
end


function value = get_field_or_empty(s, name)
    if isfield(s, name)
        value = s.(name);
    else
        value = [];
    end
end

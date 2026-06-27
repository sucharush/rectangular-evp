function spec = polygon_experiment_specs(name)
%POLYGON_EXPERIMENT_SPECS  Frozen settings for the Section 4.1 polygon experiments.
%   spec = polygon_experiment_specs(name) returns the named experiment spec.
%   polygon_experiment_specs() with no argument lists the available names.
%
%   Conventions:
%   - The QR policy is a global invariant applied by the driver, not stored per
%     spec: scan objective = CPQR with qr_tau = 1e-13; refine objective = plain
%     QR with qr_tau = []. Likewise nI = 50 and the local-refine options are
%     shared driver defaults.
%   - case_fun is a zero-argument handle returning the geometry cfg.
%   - cluster_enable = true only for the thin H-shape, which gets the
%     cluster-aware post-processing; the others are pure scan-refinement.

    specs = build_registry();

    if nargin < 1 || isempty(name)
        fprintf('Available polygon experiment specs:\n');
        fn = fieldnames(specs);
        for i = 1:numel(fn)
            fprintf('  %s\n', fn{i});
        end
        spec = specs;
        return;
    end

    if ~isfield(specs, name)
        error('polygon_experiment_specs: unknown spec ''%s''.', name);
    end
    spec = specs.(name);
end


function specs = build_registry()
    specs = struct();

    % --- pure scan-refinement experiments ---
    specs.lshape_scan = make_spec( ...
        'lshape_scan', @case_polygon_lshape, 60, [6, 20], 0.025, false);

    specs.gww1_scan = make_spec( ...
        'gww1_scan', @() case_polygon_drum('left'), 140, [1, 6], 0.025, false);

    specs.gww2_scan = make_spec( ...
        'gww2_scan', @() case_polygon_drum('right'), 140, [1, 6], 0.025, false);

    specs.hshape_width_0p3_scan = make_spec( ...
        'hshape_width_0p3_scan', @() case_polygon_hshape(0.3), 60, [0, 15], 0.025, false);

    % --- scan + cluster-aware post-processing ---
    specs.hshape_width_0p08_cluster = make_spec( ...
        'hshape_width_0p08_cluster', @() case_polygon_hshape(0.08), 60, [10, 13], 0.025, true);

    % --- GWW1-vs-GWW2 comparison: run both members and overlay their scans ---
    pair = struct();
    pair.name = 'gww_pair_scan';
    pair.pair_members = {'gww1_scan', 'gww2_scan'};
    pair.pair_labels = {'GWW1', 'GWW2'};
    pair.plot_name = 'polygon_drum_left_right_scan';
    specs.gww_pair_scan = pair;
end


function spec = make_spec(name, case_fun, disc, interval, step, cluster_enable)
    spec = struct();
    spec.name = name;
    spec.case_fun = case_fun;
    spec.Mcorner = disc;
    spec.nb_per_edge = disc;
    spec.nI = 50;
    spec.interval = interval;
    spec.step = step;
    spec.cluster_enable = cluster_enable;
    spec.plot_name = name;
end

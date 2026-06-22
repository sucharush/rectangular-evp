function h = plot_polygon_geometry(geom, varargin)
% Plot polygon vertices together with boundary and interior samples.

    p = inputParser;
    addParameter(p, 'plot_title', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'save_name', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    plot_title = char(p.Results.plot_title);
    save_name = char(p.Results.save_name);

    figure;
    h = struct();
    legend_handles = gobjects(0);
    legend_labels = {};

    h.polygon = plot(geom.V(:,1), geom.V(:,2), 'k-', 'LineWidth', 1.5);
    hold on;
    legend_handles(end+1) = h.polygon;
    legend_labels{end+1} = 'polygon';

    if isfield(geom, 'B') && ~isempty(geom.B)
        h.boundary = plot(geom.B(:,1), geom.B(:,2), 'bo', 'MarkerSize', 4);
        legend_handles(end+1) = h.boundary;
        legend_labels{end+1} = 'boundary samples';
    end

    if isfield(geom, 'I') && ~isempty(geom.I)
        h.interior = plot(geom.I(:,1), geom.I(:,2), 'r.', 'MarkerSize', 10);
        legend_handles(end+1) = h.interior;
        legend_labels{end+1} = 'interior samples';
    end

    axis equal;
    grid on;

    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'best');
    end
    if ~isempty(plot_title)
        title(plot_title, 'Interpreter', 'none');
    end
    if ~isempty(save_name)
        save_plot_eps(save_name);
    end
end

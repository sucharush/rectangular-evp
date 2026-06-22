function h = plot_singular_value_trace(lam_vec, sigma_trace, varargin)
% Plot one or more singular-value traces already computed on lam_vec.

    p = inputParser;
    addParameter(p, 'labels', {}, @(x) iscell(x) || isstring(x));
    addParameter(p, 'save_name', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    labels = p.Results.labels;
    if isstring(labels)
        labels = cellstr(labels);
    end
    save_name = char(p.Results.save_name);

    if size(sigma_trace, 1) ~= numel(lam_vec)
        if size(sigma_trace, 2) == numel(lam_vec)
            sigma_trace = sigma_trace.';
        else
            error('plot_singular_value_trace: sigma_trace must have one row per lambda.');
        end
    end

    figure;
    h = semilogy(lam_vec(:), sigma_trace, 'LineWidth', 1.5);
    grid on;
    xlabel('\lambda');
    ylabel('singular value');

    if isempty(labels)
        labels = arrayfun(@(k) sprintf('\\sigma trace %d', k), ...
            1:size(sigma_trace, 2), 'UniformOutput', false);
    end
    legend(labels, 'Location', 'best');

    if ~isempty(save_name)
        save_plot_eps(save_name);
    end
end

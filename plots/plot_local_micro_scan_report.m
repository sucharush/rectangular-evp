function h = plot_local_micro_scan_report(report, varargin)
% Plot a local cluster-resolution micro-scan report.

    p = inputParser;
    addParameter(p, 'save_name', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'plot_title', 'Local micro-scan near \lambda_*', ...
        @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    save_name = char(p.Results.save_name);
    plot_title = char(p.Results.plot_title);

    if ~isfield(report, 'micro_scan') || isempty(report.micro_scan.lamvec)
        error('plot_local_micro_scan_report: report.micro_scan is empty.');
    end

    figure;
    h = struct();
    legend_handles = gobjects(0);
    legend_labels = {};

    h.scan = semilogy(report.micro_scan.lamvec, report.micro_scan.S, ...
        'k-', 'LineWidth', 1.2);
    hold on;
    legend_handles(end+1) = h.scan;
    legend_labels{end+1} = 'micro-scan';

    grid on;
    h.lambda_star = xline(report.lam_star, 'b--', 'LineWidth', 1.0);
    legend_handles(end+1) = h.lambda_star;
    legend_labels{end+1} = '\lambda_*';

    h.base = semilogy(report.lam_star, report.sigma_min, ...
        'bo', 'MarkerFaceColor', 'b');
    legend_handles(end+1) = h.base;
    legend_labels{end+1} = 'base dip';

    extra_mask = report.micro_summary.extra_mask;
    if any(extra_mask)
        extra_lam = [report.micro_summary.refined_candidates(extra_mask).refined_lambda];
        extra_sig = [report.micro_summary.refined_candidates(extra_mask).refined_sigma];
        h.extra = semilogy(extra_lam, extra_sig, 'ro', 'MarkerFaceColor', 'r');
        legend_handles(end+1) = h.extra;
        legend_labels{end+1} = 'extra dips';
    end

    xlabel('\lambda');
    ylabel('\sigma_{min}(Q_B(\lambda))');
    % title(plot_title);
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 13);

    if ~isempty(save_name)
        save_plot_eps(save_name);
    end
end

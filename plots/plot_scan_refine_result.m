function h = plot_scan_refine_result(result, varargin)
% Plot scan curve, raw dips, and refined accepted candidates.
%
% Usage:
%   plot_scan_refine_result(result)
%   plot_scan_refine_result(result, 'problem_name', problem.name, 'save_name', 'scan_refine')
%
% Optional name-value pairs:
%   'problem_name'  : '' by default
%   'save_name'     : '' by default; pass a name to save with save_plot_eps
%   'show_raw'      : false by default
%   'show_rejected' : false by default
%   'mark_lines'    : true by default

    p = inputParser;
    addParameter(p, 'problem_name', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'save_name', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'show_raw', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'show_rejected', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'mark_lines', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    problem_name = char(p.Results.problem_name);
    save_name = char(p.Results.save_name);
    show_raw = p.Results.show_raw;
    show_rejected = p.Results.show_rejected;
    mark_lines = p.Results.mark_lines;

    if ~isfield(result, 'scan')
        error('plot_scan_refine_result: result.scan is missing.');
    end

    figure;
    h = struct();
    legend_handles = gobjects(0);
    legend_labels = {};

    h.scan = semilogy(result.scan.lamvec, result.scan.S, 'k-', 'LineWidth', 1.2);
    hold on;
    legend_handles(end+1) = h.scan;
    legend_labels{end+1} = 'scan';

    grid on;
    xlabel('\lambda');
    ylabel('\sigma_{min}(Q_B(\lambda))');
    % if ~isempty(problem_name)
    %     title(sprintf('scan-refine: %s', problem_name), 'Interpreter', 'none');
    % else
    %     title('scan-refine');
    % end

    if show_raw && isfield(result.scan, 'candidate_idx')
        J = result.scan.candidate_idx;
        if ~isempty(J)
            h.raw = plot(result.scan.lamvec(J), result.scan.S(J), ...
            'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
            legend_handles(end+1) = h.raw;
            legend_labels{end+1} = 'raw dips';
        end
    end

    if isfield(result, 'candidates') && ~isempty(result.candidates)
        cand = result.candidates;

        accepted = [cand.accepted];
        cand_acc = cand(accepted);

        if ~isempty(cand_acc)
            x_acc = [cand_acc.refined_lambda];
            y_acc = [cand_acc.refined_sigma];

            h.accepted = semilogy(x_acc, y_acc, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
            legend_handles(end+1) = h.accepted;
            legend_labels{end+1} = 'refined';

            if mark_lines
                yl = ylim;
                h.accepted_lines = gobjects(numel(cand_acc), 1);
                for i = 1:numel(cand_acc)
                    xx = cand_acc(i).refined_lambda;
                    h.accepted_lines(i) = plot([xx xx], yl, 'b--');
                end
            end
        end

        if show_rejected
            rejected = ~[cand.accepted];
            cand_rej = cand(rejected);

            keep = false(size(cand_rej));
            for i = 1:numel(cand_rej)
                keep(i) = ~isempty(cand_rej(i).refined_lambda) && ~isempty(cand_rej(i).refined_sigma) ...
                          && isfinite(cand_rej(i).refined_lambda) && isfinite(cand_rej(i).refined_sigma);
            end
            cand_rej = cand_rej(keep);

            if ~isempty(cand_rej)
                x_rej = [cand_rej.refined_lambda];
                y_rej = [cand_rej.refined_sigma];
                h.rejected = semilogy(x_rej, y_rej, 'bs', 'MarkerSize', 5);
                legend_handles(end+1) = h.rejected;
                legend_labels{end+1} = 'rejected refined';
            end
        end
    end
    % xlim([0.9, 15.1]);

    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'FontSize', 13, 'Location', 'best');
    end

    if ~isempty(save_name)
        save_plot_eps(save_name);
    end
end

function plot_scan_refine_result(result, varargin)
% Plot scan curve, raw dips, and refined accepted candidates.
%
% Usage:
%   plot_scan_refine_result(result)
%   plot_scan_refine_result(result, 'show_rejected', true)
%
% Optional name-value pairs:
%   'show_rejected' : false by default

    p = inputParser;
    addParameter(p, 'show_rejected', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    show_rejected = p.Results.show_rejected;

    if ~isfield(result, 'scan')
        error('plot_scan_refine_result: result.scan is missing.');
    end

    figure;
    plot(result.scan.lamvec, result.scan.S, 'k-', 'LineWidth', 1.2); hold on;
    grid on;
    xlabel('\lambda');
    ylabel('\sigma(\lambda)');
    title(sprintf('scan-refine: %s', result.problem_name), 'Interpreter', 'none');

    % raw scan minima
    J = result.scan.candidate_idx;
    if ~isempty(J)
        plot(result.scan.lamvec(J), result.scan.S(J), ...
            'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    end

    if isfield(result, 'candidates') && ~isempty(result.candidates)
        cand = result.candidates;

        accepted = [cand.accepted];
        cand_acc = cand(accepted);

        if ~isempty(cand_acc)
            x_acc = [cand_acc.refined_lambda];
            y_acc = [cand_acc.refined_sigma];

            plot(x_acc, y_acc, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

            yl = ylim;
            for i = 1:numel(cand_acc)
                xx = cand_acc(i).refined_lambda;
                plot([xx xx], yl, 'r--');
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
                plot(x_rej, y_rej, 'bs', 'MarkerSize', 5);
            end
        end
    end

    lgd = {'scan', 'raw dips', 'accepted refined'};
    if show_rejected
        lgd{end+1} = 'rejected refined';
    end
    legend(lgd, 'Location', 'best');
end
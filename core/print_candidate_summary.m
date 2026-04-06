function print_candidate_summary(result)
% Pretty print candidate records from solve_scan_refine output.

    if ~isfield(result, 'candidates') || isempty(result.candidates)
        fprintf('No candidates available.\n');
        return;
    end

    cand = result.candidates;

    % fprintf('\n=== Candidate summary: %s ===\n', result.problem_name);
    fprintf('\n=== Candidate summary: ===\n');
    fprintf('%4s  %8s  %14s  %14s  %14s  %10s  %8s\n', ...
        'id', 'status', 'scan_lambda', 'refined_lambda', 'refined_sigma', 'n_calls', 'bracket');

    for k = 1:numel(cand)
        ck = cand(k);

        id_str = sprintf('%d', ck.id);

        if isempty(ck.status)
            status_str = '-';
        else
            status_str = ck.status;
        end

        if isempty(ck.scan_lambda)
            scan_lam_str = '-';
        else
            scan_lam_str = sprintf('%.8f', ck.scan_lambda);
        end

        if isempty(ck.refined_lambda)
            ref_lam_str = '-';
        else
            ref_lam_str = sprintf('%.12f', ck.refined_lambda);
        end

        if isempty(ck.refined_sigma)
            ref_sig_str = '-';
        else
            ref_sig_str = sprintf('%.3e', ck.refined_sigma);
        end

        if isfield(ck, 'logs') && isfield(ck.logs, 'n_calls') && ~isempty(ck.logs.n_calls)
            ncall_str = sprintf('%d', ck.logs.n_calls);
        else
            ncall_str = '-';
        end

        if isempty(ck.bracket)
            bracket_str = '-';
        else
            bracket_str = sprintf('[%.4f,%.4f]', ck.bracket(1), ck.bracket(2));
        end

        fprintf('%4s  %8s  %14s  %14s  %14s  %10s  %8s\n', ...
            id_str, status_str, scan_lam_str, ref_lam_str, ref_sig_str, ncall_str, bracket_str);

        if isfield(ck, 'logs') && isfield(ck.logs, 'message') && ~isempty(ck.logs.message)
            fprintf('      note: %s\n', ck.logs.message);
        end
    end

    fprintf('\nAccepted: %d / %d\n', result.summary.n_accepted, result.summary.n_total);
end
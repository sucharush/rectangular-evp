function summary = summarize_candidates(candidates)
% Build a compact summary from candidate records.
%
% Keeps the detailed candidate structs intact elsewhere.

    if isempty(candidates)
        summary = struct();
        summary.n_total = 0;
        summary.n_accepted = 0;
        summary.accepted_mask = false(0,1);
        summary.eigs = [];
        summary.sigmins = [];
        summary.brackets = zeros(0,2);
        summary.scan_indices = [];
        summary.status_list = {};
        return;
    end

    accepted_mask = [candidates.accepted].';
    idx = find(accepted_mask);

    summary = struct();
    summary.n_total = numel(candidates);
    summary.n_accepted = numel(idx);
    summary.accepted_mask = accepted_mask;
    summary.status_list = {candidates.status}.';

    if isempty(idx)
        summary.eigs = [];
        summary.sigmins = [];
        summary.brackets = zeros(0,2);
        summary.scan_indices = [];
    else
        summary.eigs = [candidates(idx).refined_lambda].';
        summary.sigmins = [candidates(idx).refined_sigma].';
        summary.brackets = vertcat(candidates(idx).bracket);
        summary.scan_indices = [candidates(idx).scan_index].';
    end
end
% ==============================================================================
% 04a_extract_edges.m
%
% Extract significant edges (with t-values) from NBS result .mat files as CSV.
% Walks one or more directory trees, finds every nbs_result_*.mat, and writes
% nbs_result_*_significant_edges.csv next to each source file.
%
% Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
%           in atrial fibrillation" (Akdeniz et al.)
% Methods: Network-based statistics (post-processing of NBS results;
%          Supplementary Table 8 lists the small Control > AF component)
%
% Output CSV columns:
%   Edge_ID, Region1, Region2, t_value, abs_t
%
% Edge counts use the upper triangle of the symmetric component mask
% (no double-counting). t-values are read from nbs.NBS.test_stat with
% fallback to nbs.STATS.test_stat for older NBS versions.
% ==============================================================================
clear; clc;

% ------------------------------------------------------------------------------
% Configuration  --- edit paths for your environment
% ------------------------------------------------------------------------------
RESULTS_ROOTS = {
    '/path/to/nbs/02_results/fc'   % FC NBS results
    '/path/to/nbs/02_results/sc'   % SC NBS results
};

% Optional filters (leave {} to include everything)
include_keywords = {};   % e.g. {'design_1_total','contrast_01_af_gt_control','T4.5'}
exclude_keywords = {};

% ------------------------------------------------------------------------------
% Discover all nbs_result_*.mat under the result roots
% ------------------------------------------------------------------------------
mat_list = {};

for r = 1:numel(RESULTS_ROOTS)
    rootDir = RESULTS_ROOTS{r};
    if ~isfolder(rootDir)
        warning('Root not found (skipping): %s', rootDir);
        continue;
    end

    hits = dir(fullfile(rootDir, '**', 'nbs_result_*.mat'));   % R2016b+
    for k = 1:numel(hits)
        fp = fullfile(hits(k).folder, hits(k).name);
        if all_keywords_match(fp, include_keywords) && ...
           ~any_keyword_matches(fp, exclude_keywords)
            mat_list{end + 1, 1} = fp; %#ok<SAGROW>
        end
    end
end

mat_list = unique(mat_list);
fprintf('Discovered %d NBS result files.\n\n', numel(mat_list));

if isempty(mat_list)
    fprintf('Nothing to do.\n');
    return;
end

% ------------------------------------------------------------------------------
% Process each file and export edge CSV
% ------------------------------------------------------------------------------
fprintf('=== Batch NBS edge export (%d files) ===\n\n', numel(mat_list));

for f = 1:numel(mat_list)
    nbs_file_path = mat_list{f};
    fprintf('(%d/%d) %s\n', f, numel(mat_list), nbs_file_path);

    S = load(nbs_file_path, 'nbs');
    if ~isfield(S, 'nbs')
        warning('  No ''nbs'' variable. Skipping.');
        continue;
    end
    nbs = S.nbs; clear S;

    try
        % --- FWER p-value and permutation count ---
        if isfield(nbs, 'NBS') && isfield(nbs.NBS, 'pval') && ~isempty(nbs.NBS.pval)
            p_value = nbs.NBS.pval(1);
        else
            p_value = NaN;
        end
        if isfield(nbs, 'GLM') && isfield(nbs.GLM, 'perms')
            num_perms = nbs.GLM.perms;
        else
            num_perms = NaN;
        end

        % --- Binary component mask ---
        if isfield(nbs, 'NBS') && isfield(nbs.NBS, 'con_mat') && ~isempty(nbs.NBS.con_mat)
            network_mat = full(nbs.NBS.con_mat{1}) > 0;
        else
            error('Missing NBS.con_mat{1}');
        end

        % --- Node labels ---
        if isfield(nbs, 'NBS') && isfield(nbs.NBS, 'node_label') && ~isempty(nbs.NBS.node_label)
            node_names = nbs.NBS.node_label;
        elseif isfield(nbs, 'node_label') && ~isempty(nbs.node_label)
            node_names = nbs.node_label;
        else
            error('Missing node labels');
        end
        if ischar(node_names), node_names = cellstr(node_names); end

        % --- Test-stat matrix (t-values) ---
        if isfield(nbs, 'NBS') && isfield(nbs.NBS, 'test_stat') && ~isempty(nbs.NBS.test_stat)
            Tmat = full(nbs.NBS.test_stat);
        elseif isfield(nbs, 'STATS') && isfield(nbs.STATS, 'test_stat') && ~isempty(nbs.STATS.test_stat)
            Tmat = full(nbs.STATS.test_stat);
        else
            error('Could not find test-stat matrix');
        end

        if ~isequal(size(network_mat), size(Tmat))
            error('Component mask and t-stat matrix size mismatch');
        end

    catch ME
        warning('  Unexpected NBS structure (%s). Skipping.', ME.message);
        continue;
    end

    % --- Edge enumeration (upper triangle) ---
    [i1, i2] = find(triu(network_mat, 1));
    num_edges = numel(i1);

    if ~isnan(p_value)
        if p_value == 0 && ~isnan(num_perms) && num_perms > 0
            fprintf('  p < %g (FWER, %d perms), %d edges\n', 1/max(1, num_perms), num_perms, num_edges);
        else
            fprintf('  p = %g (FWER), %d edges\n', p_value, num_edges);
        end
    else
        fprintf('  p not available, %d edges\n', num_edges);
    end

    % --- Write CSV next to the .mat ---
    [fp, nm, ~] = fileparts(nbs_file_path);
    out_csv = fullfile(fp, [nm, '_significant_edges.csv']);

    if num_edges > 0
        lin    = sub2ind(size(Tmat), i1, i2);
        t_vals = Tmat(lin);
        edges_table = table( ...
            (1:num_edges).', node_names(i1), node_names(i2), t_vals, abs(t_vals), ...
            'VariableNames', {'Edge_ID','Region1','Region2','t_value','abs_t'});
    else
        edges_table = table([], {}, {}, [], [], ...
            'VariableNames', {'Edge_ID','Region1','Region2','t_value','abs_t'});
    end

    try
        writetable(edges_table, out_csv);
        fprintf('  -> %s\n\n', out_csv);
    catch ME
        warning('  CSV write failed: %s\n', ME.message);
    end
end

fprintf('=== Done. ===\n');

% ------------------------------------------------------------------------------
% Local helpers
% ------------------------------------------------------------------------------
function ok = all_keywords_match(s, kws)
    ok = true;
    for k = 1:numel(kws)
        if ~contains(s, kws{k}), ok = false; return; end
    end
end

function ok = any_keyword_matches(s, kws)
    ok = false;
    for k = 1:numel(kws)
        if contains(s, kws{k}), ok = true; return; end
    end
end

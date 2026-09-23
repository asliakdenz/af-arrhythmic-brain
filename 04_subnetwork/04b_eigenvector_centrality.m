% ==============================================================================
% 04b_eigenvector_centrality.m
%
% Compute eigenvector centrality (EC) on each significant NBS subnetwork,
% in both binary (presence/absence of edge) and weighted (|t-value|) form.
% Walks one or more directory trees, finds every nbs_result_*.mat, and writes
% two CSVs next to each source file:
%
%   <basename>_binary_hubs_eigen.csv     RegionName, BinaryEigenvectorCentrality
%   <basename>_weighted_hubs_eigen.csv   RegionName, WeightedEigenvectorCentrality
%
% Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
%           in atrial fibrillation" (Akdeniz et al.)
% Methods: Network-based statistics -> Functional NBS and hub analysis
%          (weighted eigenvector centrality with |t| edge weights; Figure 2;
%          Supplementary Table 9 after 04e_harvard_oxford_labels.py)
%
% Dependencies:
%   - Brain Connectivity Toolbox (Rubinov & Sporns 2010), eigenvector_centrality_und
%   - BrainNet Viewer is used separately to render Figure 2 from the hub CSV
%
% Inputs:
%   - NBS result .mat files under RESULTS_ROOTS
%   - 216-node atlas labels file (canonical node ordering)
% ==============================================================================
clear; clc;

% ------------------------------------------------------------------------------
% Configuration  --- edit paths for your environment
% ------------------------------------------------------------------------------
BCT_PATH   = '/path/to/BCT';

LABEL_FILE = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'atlas', 'labels_216.txt');

RESULTS_ROOTS = {
    '/path/to/nbs/02_results/fc'
    '/path/to/nbs/02_results/sc'
};

include_keywords = {};   % e.g. {'design_1_total','contrast_01_af_gt_control'}
exclude_keywords = {};

% ------------------------------------------------------------------------------
% Setup
% ------------------------------------------------------------------------------
addpath(genpath(BCT_PATH));

fprintf('--- Eigenvector hub analysis (binary + weighted) ---\n');
node_names_override = read_label_list(LABEL_FILE);
fprintf('Atlas: %d nodes from %s\n\n', numel(node_names_override), LABEL_FILE);

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
    hits = dir(fullfile(rootDir, '**', 'nbs_result_*.mat'));
    for k = 1:numel(hits)
        fp = fullfile(hits(k).folder, hits(k).name);
        if all_keywords_match(fp, include_keywords) && ...
           ~any_keyword_matches(fp, exclude_keywords)
            mat_list{end + 1, 1} = fp; %#ok<SAGROW>
        end
    end
end
mat_list = unique(mat_list);
fprintf('Discovered %d NBS files.\n\n', numel(mat_list));
if isempty(mat_list), return; end

% ------------------------------------------------------------------------------
% Process each .mat and save hub CSVs next to it
% ------------------------------------------------------------------------------
for f = 1:numel(mat_list)
    fp         = mat_list{f};
    folder_dir = fileparts(fp);
    fprintf('(%d/%d) %s\n', f, numel(mat_list), fp);

    S = load(fp, 'nbs');
    if ~isfield(S, 'nbs')
        warning('  No ''nbs'' variable. Skipping.'); continue;
    end
    nbs = S.nbs; clear S;

    % --- Component mask ---
    if ~isfield(nbs, 'NBS') || ~isfield(nbs.NBS, 'con_mat') || isempty(nbs.NBS.con_mat)
        fprintf('  No significant component. Skipping.\n'); continue;
    end
    cm  = nbs.NBS.con_mat;
    Aup = iscell(cm) * 0;   % placeholder
    if iscell(cm)
        if isempty(cm{1})
            fprintf('  Empty component. Skipping.\n'); continue;
        end
        Aup = cm{1};
    else
        Aup = cm;
    end
    if isempty(Aup) || nnz(Aup) == 0
        fprintf('  Zero edges. Skipping.\n'); continue;
    end

    % --- t-statistics for weighted EC ---
    T = [];
    if isfield(nbs, 'NBS') && isfield(nbs.NBS, 'test_stat') && ~isempty(nbs.NBS.test_stat)
        ts = nbs.NBS.test_stat;
        if iscell(ts), ts = ts{1}; end
        T = full(ts);
    elseif isfield(nbs, 'STATS') && isfield(nbs.STATS, 'test_stat') && ~isempty(nbs.STATS.test_stat)
        ts = nbs.STATS.test_stat;
        if iscell(ts), ts = ts{1}; end
        T = full(ts);
    end

    % --- Dimension check ---
    N = size(Aup, 1);
    if numel(node_names_override) ~= N
        warning('  Label count (%d) != adjacency size (%d). Skipping.', ...
                numel(node_names_override), N);
        continue;
    end
    node_names = node_names_override;

    % --- Binary symmetric adjacency ---
    Abin = double(Aup > 0);
    Abin = Abin + Abin.';
    Abin(1:N + 1:end) = 0;

    % --- Weighted symmetric adjacency (|t| on present edges) ---
    if ~isempty(T)
        W = abs(T) .* double(Aup > 0);
    else
        W = double(Aup > 0);
        warning('  test_stat missing; weighted EC reduces to binary EC.');
    end
    W = W + W.';
    W(1:N + 1:end) = 0;

    % --- Eigenvector centrality (BCT) ---
    ec_bin = eigenvector_centrality_und(full(Abin));
    ec_wtd = eigenvector_centrality_und(full(W));

    % --- Sort and filter EC > 0 ---
    [ec_bin_sorted, idx] = sort(ec_bin, 'descend');
    names_sorted  = node_names(idx);
    ec_wtd_sorted = ec_wtd(idx);
    keep_bin = ec_bin_sorted > 0;
    keep_wtd = ec_wtd_sorted > 0;

    T_bin = table(names_sorted(keep_bin), ec_bin_sorted(keep_bin), ...
        'VariableNames', {'RegionName','BinaryEigenvectorCentrality'});
    T_wtd = table(names_sorted(keep_wtd), ec_wtd_sorted(keep_wtd), ...
        'VariableNames', {'RegionName','WeightedEigenvectorCentrality'});

    [~, base, ~] = fileparts(fp);
    out_bin = fullfile(folder_dir, [base, '_binary_hubs_eigen.csv']);
    out_wtd = fullfile(folder_dir, [base, '_weighted_hubs_eigen.csv']);

    try
        writetable(T_bin, out_bin);
        writetable(T_wtd, out_wtd);
        fprintf('  -> %d binary hubs, %d weighted hubs (EC > 0)\n\n', ...
                height(T_bin), height(T_wtd));
    catch ME
        warning('  Write failed: %s', ME.message);
    end
end

fprintf('--- Done. Outputs saved next to each nbs_result_*.mat. ---\n');

% ------------------------------------------------------------------------------
% Local helpers
% ------------------------------------------------------------------------------
function labels = read_label_list(txtfile)
    if ~isfile(txtfile)
        error('Label file not found: %s', txtfile);
    end
    raw = readlines(txtfile);
    raw = strtrim(raw);
    raw(raw == "") = [];
    labels = cellstr(raw);
    labels = labels(:);
end

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

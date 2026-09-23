% ==============================================================================
% 04c_threshold_summary.m
%
% Stability of an NBS solution across primary t-thresholds: component size
% and FWER p-value per threshold (Supplementary Tables 5a and 6), and overlap
% of the higher-threshold components with the primary (lowest-threshold)
% solution (Supplementary Table 5b).
%
% Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
%           in atrial fibrillation" (Akdeniz et al.)
% Methods: Network-based statistics -> Functional NBS and hub analysis
%          ("We repeated the analysis from t = 4.5 to 6.0 in steps of 0.5,
%          reporting component significance and extent at each threshold")
%
% Overlap measures (Supplementary Table 5b), each threshold versus the
% primary solution:
%   edge Jaccard      |A n B| / |A u B|
%   edge containment  |A n B| / min(|A|, |B|)   (1 = strictly nested)
%   node Jaccard      on the node sets of the components
%   hub-rank rho      Spearman correlation of weighted eigenvector-centrality
%                     ranks across the nodes shared by the two solutions
%                     (weighted hub CSVs from 04b_eigenvector_centrality.m),
%                     with a t-approximation P-value
%
% Expected folder structure (output of 03d/03e, 04a and 04b):
%   <baseDir>/T4.5/Extent/nbs_result_T4.5_Extent.mat
%   <baseDir>/T4.5/Extent/nbs_result_T4.5_Extent_significant_edges.csv
%   <baseDir>/T4.5/Extent/nbs_result_T4.5_Extent_weighted_hubs_eigen.csv
%   <baseDir>/T5/Extent/ ...
%
% Outputs:
%   <baseDir>/nbs_threshold_summary_extent.csv
%   <baseDir>/nbs_threshold_stability.csv
%
% Run once per (design x contrast) directory.
% ==============================================================================
clear; clc;

% ------------------------------------------------------------------------------
% Configuration  --- set the directory containing the T*/Extent subfolders
% ------------------------------------------------------------------------------
baseDir = '/path/to/nbs/02_results/fc/design_1_total/contrast_01_af_gt_control';

outCsv       = fullfile(baseDir, 'nbs_threshold_summary_extent.csv');
stabilityCsv = fullfile(baseDir, 'nbs_threshold_stability.csv');

% ------------------------------------------------------------------------------
% Thresholds from T* subfolders
% ------------------------------------------------------------------------------
d = dir(fullfile(baseDir, 'T*')); d = d([d.isdir]);
tVals = [];
for i = 1:numel(d)
    tok = regexp(d(i).name, '^T([0-9]+(\.[0-9]+)?)$', 'tokens', 'once');
    if ~isempty(tok), tVals(end + 1) = str2double(tok{1}); end %#ok<SAGROW>
end
if isempty(tVals), error('No threshold folders found under: %s', baseDir); end
tVals = sort(unique(tVals));

rows = struct('t', {}, 'Nodes', {}, 'Edges', {}, 'p_FWER', {}, 'p_fmt', {}, ...
              'csv_path', {}, 'mat_path', {});
edgeSets = cell(numel(tVals), 1);   % sorted "A|B" edge keys per threshold
nodeSets = cell(numel(tVals), 1);
hubTabs  = cell(numel(tVals), 1);   % table(RegionName, EC) per threshold

% ------------------------------------------------------------------------------
% Per-threshold summary
% ------------------------------------------------------------------------------
for k = 1:numel(tVals)
    t    = tVals(k);
    tStr = regexprep(num2str(t), '\.0$', '');
    extentDir = fullfile(baseDir, ['T' tStr], 'Extent');

    csvPath = first_match(extentDir, sprintf('nbs_result_T%s_Extent_significant_edges.csv', tStr), '*significant_edges*.csv');
    matPath = first_match(extentDir, sprintf('nbs_result_T%s_Extent.mat', tStr), '*.mat');
    hubPath = first_match(extentDir, sprintf('nbs_result_T%s_Extent_weighted_hubs_eigen.csv', tStr), '*weighted_hubs_eigen.csv');

    nEdges = 0; nNodes = 0; keys = strings(0, 1); nodes = strings(0, 1);
    if ~isempty(csvPath)
        T = readtable(csvPath, 'TextType', 'string');
        nEdges = height(T);
        if nEdges > 0
            r1 = string(T.Region1); r2 = string(T.Region2);
            keys = strings(nEdges, 1);
            for e = 1:nEdges
                pair = sort([r1(e) r2(e)]);
                keys(e) = pair(1) + "|" + pair(2);
            end
            nodes  = unique([r1; r2]);
            nNodes = numel(nodes);
        end
    end
    edgeSets{k} = keys; nodeSets{k} = nodes;

    if ~isempty(hubPath)
        H = readtable(hubPath, 'TextType', 'string');
        ecCol = H.Properties.VariableNames(contains(H.Properties.VariableNames, 'Eigenvector'));
        if ~isempty(ecCol), hubTabs{k} = table(string(H.RegionName), H.(ecCol{1}), 'VariableNames', {'RegionName', 'EC'}); end
    end

    pFWER = NaN; pLowerBound = NaN; numPerms = NaN;
    if ~isempty(matPath)
        S = load(matPath);
        if isfield(S, 'nbs') && isstruct(S.nbs)
            if isfield(S.nbs, 'GLM') && isfield(S.nbs.GLM, 'perms'), numPerms = double(S.nbs.GLM.perms); end
            if isfield(S.nbs, 'NBS') && isfield(S.nbs.NBS, 'pval'), pFWER = tryScalarP(S.nbs.NBS.pval); end
        end
        % p == 0 from a permutation test means p < 1 / permutations
        if ~isnan(pFWER) && pFWER == 0 && ~isnan(numPerms) && numPerms > 0, pLowerBound = 1 / numPerms; end
    end

    rows(k).t = t;
    rows(k).Nodes = string(nNodes); rows(k).Edges = string(nEdges);
    if nEdges == 0, rows(k).Nodes = "-"; rows(k).Edges = "-"; end
    rows(k).p_FWER = pFWER; rows(k).p_fmt = formatP(pFWER, pLowerBound);
    rows(k).csv_path = string(csvPath); rows(k).mat_path = string(matPath);
end

Tout = table([rows.t]', string({rows.Nodes})', string({rows.Edges})', [rows.p_FWER]', ...
             string({rows.p_fmt})', string({rows.csv_path})', string({rows.mat_path})', ...
    'VariableNames', {'Primary_threshold_t', 'Nodes', 'Edges', 'p_FWER', 'p_FWER_formatted', 'csv_path', 'mat_path'});
writetable(Tout, outCsv);
disp(['Saved: ' outCsv]);
disp(Tout(:, {'Primary_threshold_t', 'Nodes', 'Edges', 'p_FWER_formatted'}));

% ------------------------------------------------------------------------------
% Stability versus the primary (lowest) threshold  (Supplementary Table 5b)
% ------------------------------------------------------------------------------
if numel(tVals) > 1 && ~isempty(edgeSets{1})
    A = edgeSets{1}; nodesA = nodeSets{1};
    stab = table('Size', [numel(tVals) - 1, 7], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Comparison', 'Edge_Jaccard', 'Edge_containment', 'Node_Jaccard', ...
                          'Hub_rank_rho', 'n_shared_nodes', 'Hub_rank_P'});
    for k = 2:numel(tVals)
        B = edgeSets{k}; nodesB = nodeSets{k};
        stab.Comparison(k - 1) = sprintf('t = %g vs %g', tVals(1), tVals(k));
        if isempty(B), continue; end
        inter = numel(intersect(A, B));
        stab.Edge_Jaccard(k - 1)     = inter / numel(union(A, B));
        stab.Edge_containment(k - 1) = inter / min(numel(A), numel(B));
        stab.Node_Jaccard(k - 1)     = numel(intersect(nodesA, nodesB)) / numel(union(nodesA, nodesB));
        if ~isempty(hubTabs{1}) && ~isempty(hubTabs{k})
            [shared, ia, ib] = intersect(hubTabs{1}.RegionName, hubTabs{k}.RegionName);
            n = numel(shared);
            stab.n_shared_nodes(k - 1) = n;
            if n >= 3
                [rho, p] = spearman_local(hubTabs{1}.EC(ia), hubTabs{k}.EC(ib));
                stab.Hub_rank_rho(k - 1) = rho; stab.Hub_rank_P(k - 1) = p;
            end
        end
    end
    writetable(stab, stabilityCsv);
    disp(['Saved: ' stabilityCsv]);
    disp(stab);
end

% ==============================================================================
% Local helpers
% ==============================================================================
function p = first_match(folder, preferred, pattern)
    p = '';
    if isfile(fullfile(folder, preferred)), p = fullfile(folder, preferred); return; end
    c = dir(fullfile(folder, pattern));
    if ~isempty(c), p = fullfile(c(1).folder, c(1).name); end
end

function p = tryScalarP(x)
    p = NaN;
    try
        x = double(x); x = x(:);
        if ~isempty(x) && isfinite(x(1)) && x(1) >= 0 && x(1) <= 1, p = x(1); end
    catch
    end
end

function s = formatP(p, pLowerBound)
    if isnan(p), s = "NA"; return; end
    if p == 0 && ~isnan(pLowerBound)
        s = sprintf('< %.4g', pLowerBound);
    else
        s = sprintf('%.4f', p);
    end
end

function rk = tiedrank_local(x)
    x = x(:); [s, idx] = sort(x); n = numel(x); rk = zeros(n, 1); i = 1;
    while i <= n
        j = i;
        while j < n && s(j + 1) == s(i), j = j + 1; end
        rk(idx(i:j)) = (i + j) / 2;
        i = j + 1;
    end
end

function [rho, p] = spearman_local(x, y)
    rx = tiedrank_local(x); ry = tiedrank_local(y);
    rx = rx - mean(rx); ry = ry - mean(ry);
    rho = (rx' * ry) / sqrt((rx' * rx) * (ry' * ry));
    n = numel(x); df = n - 2;
    if abs(rho) >= 1, p = 0; return; end
    tstat = rho * sqrt(df / (1 - rho^2));
    p = betainc(df / (df + tstat^2), df / 2, 0.5);   % two-sided t-approximation
end

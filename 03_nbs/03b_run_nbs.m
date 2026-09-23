% ==============================================================================
% 03b_run_nbs.m
%
% MATLAB runner for a single NBS analysis using NBS-Connectome v1.2
% (Zalesky et al., 2010).
%
% Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
%           in atrial fibrillation" (Akdeniz et al.)
% Methods: Network-based statistics
%
% Inputs (strings, passed from 03c_run_nbs.sh):
%   output_dir     directory where nbs_result_T*.mat will be written
%   matrices_path  path to any matrix file in the matrix directory
%                  (NBS reads every file of the directory, sorted alphabetically)
%   design_path    space-delimited design.txt (no header)
%   contrast_path  space-delimited single-row contrast vector
%   thresh_str     primary t-threshold as a string (e.g. '4.5')
%   size_str       component size statistic: 'Extent' (used) or 'Intensity'
%   method_str     NBS method ('Run NBS' for standard NBS)
%
% Fixed settings (Methods): 5,000 permutations, alpha = 0.05, t-test,
% 216-node Schaefer 200 + Tian S1 atlas (data/atlas of this repository).
% ==============================================================================
function run_nbs(output_dir, matrices_path, design_path, contrast_path, ...
                 thresh_str, size_str, method_str)

    fprintf('--- Starting NBS Analysis ---\n');

    % NBS toolbox path (edit for your environment)
    NBS_TOOLBOX_PATH = '/path/to/NBS1.2/';
    addpath(genpath(NBS_TOOLBOX_PATH));

    fprintf('Output Dir:  %s\n', output_dir);
    fprintf('Matrix File: %s\n', matrices_path);
    fprintf('Design:      %s\n', design_path);
    fprintf('Contrast:    %s\n', contrast_path);
    fprintf('Method:      %s\n', method_str);
    fprintf('Threshold:   %s\n', thresh_str);
    fprintf('Size:        %s\n', size_str);

    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    try, maxNumCompThreads(2); catch, end
    clear global nbs; global nbs;

    % --- Atlas files: data/atlas of this repository ---
    here      = fileparts(mfilename('fullpath'));
    atlas_dir = fullfile(here, '..', 'data', 'atlas');
    node_coords_path = fullfile(atlas_dir, 'coords_216.txt');
    node_labels_path = fullfile(atlas_dir, 'labels_216.txt');
    if ~exist(node_coords_path, 'file'), error('Missing node coordinates file: %s', node_coords_path); end
    if ~exist(node_labels_path, 'file'), error('Missing node labels file: %s', node_labels_path); end

    % --- NBS UI structure ---
    UI.test.ui       = 't-test';
    UI.perms.ui      = '5000';
    UI.alpha.ui      = '0.05';
    UI.exchange.ui   = '';
    UI.method.ui     = method_str;
    UI.size.ui       = canonical_size(size_str);
    UI.thresh.ui     = thresh_str;
    UI.contrast.ui   = contrast_path;
    UI.design.ui     = design_path;
    UI.matrices.ui   = matrices_path;
    UI.node_coor.ui  = node_coords_path;
    UI.node_label.ui = node_labels_path;

    disp('Calling NBSrun...');
    NBSrun(UI, []);

    if isfield(nbs, 'NBS') && isfield(nbs.NBS, 'n') && nbs.NBS.n > 0
        outname = sprintf('nbs_result_T%s_%s.mat', UI.thresh.ui, UI.size.ui);
        save(fullfile(output_dir, outname), 'nbs');
        fprintf('--- Significant networks found and saved to %s ---\n', outname);
    else
        fprintf('--- No significant networks found.\n');
    end
    fprintf('--- Finished ---\n');
    exit;
end

function s = canonical_size(x)
    if strcmpi(x, 'extent'),    s = 'Extent';    return; end
    if strcmpi(x, 'intensity'), s = 'Intensity'; return; end
    s = 'Extent';
end

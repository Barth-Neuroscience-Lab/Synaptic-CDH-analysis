%% ============================================================
%  BALANCED SUBSAMPLING + SUBJECT-LEVEL PERMUTATION
%  TWO-SAMPLE KOLMOGOROV-SMIRNOV (K-S) ANALYSIS
%  ============================================================
%  File: balanced_subsample_perm_ks_analysis_v2.m
%
%  VERSION HISTORY
%    v1  Balanced within-subject subsampling, subject-level permutation
%        K-S test, uncorrected and Benjamini-Hochberg p-values.
%    v2  Adds Holm-Bonferroni and Westfall-Young (min-P, global subject
%        permutation) family-wise corrections; PLOT_CORRECTION option.
%
%  PURPOSE
%    Compare the distribution of a measurement between experimental
%    groups when each group contains several subjects and each subject
%    contributes many observations. Every subject is given equal weight
%    by drawing the same number of observations from each one.
%
%  INPUT
%    An Excel file chosen via a dialog. Each sheet = one group.
%    Each column = one subject (name in the first non-empty cell,
%    numeric observations below). Non-numeric / non-finite cells are
%    ignored.
%
%  PROCEDURE
%    1. Balanced subsampling (observed data, real group labels)
%       For each of N_ITER iterations:
%         - From EVERY subject in EVERY selected group, draw N_SAMP
%           observations WITHOUT replacement. N_SAMP is identical for
%           all subjects and cannot exceed the smallest per-subject
%           sample count among the selected groups.
%         - Pool the draws within each group.
%         - For every pair of groups, compute the two-sample K-S
%           D statistic, D = max |F1(x) - F2(x)|.
%         - Record each group's empirical CDF on a fixed grid.
%       -> Distribution of observed D for each pair (N_ITER values).
%
%    2. Subject-level permutation (null distribution)
%       For each pair and each of N_PERM permutations:
%         - Pool the subjects of both groups and randomly reassign
%           them (without replacement) to two pseudo-groups with the
%           original numbers of subjects.
%         - Draw N_SAMP observations per subject exactly as in step 1.
%         - Compute D.
%       -> Null distribution of D for each pair (N_PERM values).
%
%    3. Inference
%       Uncorrected p = (1 + #{null D >= mean observed D}) / (1 + N_PERM)
%       Three multiple-comparison corrections across pairs:
%         - Benjamini-Hochberg (controls the false discovery rate)
%         - Holm-Bonferroni   (controls the family-wise error rate)
%         - Westfall-Young single-step max-statistic (min-P) permutation
%           correction (controls the family-wise error rate; the
%           nonparametric counterpart of Tukey's HSD). Subjects are
%           permuted across ALL selected groups at once (group sizes
%           kept), N_SAMP observations are drawn per subject, and D is
%           computed for every pair. Each pair's D is converted to its
%           null p-value so pairs with different numbers of subjects are
%           on a common scale, and the smallest p across pairs is kept
%           for each permutation. A pair's corrected p is the fraction
%           of permutations whose smallest p is <= that pair's p.
%
%  OUTPUTS (saved next to the input file)
%    <name>_BalSubPerm_KS_results.xlsx with sheets:
%      Parameters, Subject_Statistics, KS_Summary,
%      Observed_D_Iterations, Permutation_D_Null,
%      Global_Permutation_D_Null (3+ groups), Mean_CDF_by_Pair
%    PNG figures: observed-vs-null D per pair, mean CDFs per pair,
%      and a summary of mean observed D vs the null critical value.
% ============================================================

clearvars; close all; clc;

%% ── FIXED CONFIG ────────────────────────────────────────────
ALPHA = 0.05;     % significance level for all corrections
PLOT_CORRECTION = 'BH';  % correction marked with * on the summary plot:
                         % 'BH', 'Holm' or 'WY' (Westfall-Young max-D)
SEED  = 42;       % set to 0 for a non-reproducible random stream
MAX_BINS = 2000;  % upper limit on CDF grid points (memory guard)
%% ─────────────────────────────────────────────────────────────

if SEED ~= 0, rng(SEED); end

%% ── 1. USER SELECTS INPUT FILE ──────────────────────────────
[fname, fpath] = uigetfile({'*.xlsx;*.xls', 'Excel Files (*.xlsx, *.xls)'}, ...
                            'Select your data file');
if isequal(fname, 0)
    error('No file selected. Script aborted.');
end
EXCEL_FILE = fullfile(fpath, fname);

[~, baseName, ~] = fileparts(fname);
OUT_TAG     = '_BalSubPerm_KS';
OUTPUT_FILE = fullfile(fpath, [baseName OUT_TAG '_results.xlsx']);

fprintf('Input file : %s\n', EXCEL_FILE);
fprintf('Output file: %s\n\n', OUTPUT_FILE);

%% ── 2. READ DATA ─────────────────────────────────────────────
fprintf('Reading data ...\n');
sheets  = cellstr(sheetnames(EXCEL_FILE));
nGroups = numel(sheets);
fprintf('  Found %d group(s): %s\n\n', nGroups, strjoin(sheets, ', '));

groupData      = cell(nGroups, 1);
groupSubjNames = cell(nGroups, 1);

for g = 1:nGroups
    raw = readcell(EXCEL_FILE, 'Sheet', sheets{g});
    [nRows, nCols] = size(raw);
    subjects  = cell(nCols, 1);
    subjNames = cell(nCols, 1);
    nSubjects = 0;

    for c = 1:nCols
        col = raw(:, c);
        firstNonEmpty = find(~cellfun(@is_blank_cell, col), 1);
        if isempty(firstNonEmpty), continue; end

        rawName = col{firstNonEmpty};
        if ischar(rawName) || isstring(rawName)
            sName = char(strtrim(string(rawName)));
        elseif isnumeric(rawName) && isscalar(rawName)
            sName = num2str(rawName);
        else
            sName = sprintf('Subject_%d', c);
        end

        dataStart = firstNonEmpty + 1;
        if dataStart > nRows, continue; end

        vals   = col(dataStart:end);
        isNum  = cellfun(@(v) isnumeric(v) && isscalar(v), vals);
        numVec = cell2mat(vals(isNum));
        numVec = double(numVec(isfinite(numVec)));
        if isempty(numVec), continue; end

        nSubjects = nSubjects + 1;
        subjects{nSubjects}  = numVec(:);
        subjNames{nSubjects} = sName;
    end

    groupData{g}      = subjects(1:nSubjects);
    groupSubjNames{g} = subjNames(1:nSubjects);
end

%% ── 3. EXCLUDE GROUPS WITH < 2 SUBJECTS ─────────────────────
% Subject-level permutation needs at least 2 subjects per group.
nSubjAll = cellfun(@numel, groupData);
tooFew   = nSubjAll < 2;
for g = find(tooFew)'
    fprintf('  WARNING: Group "%s" has %d valid subject(s) and is excluded\n', ...
            sheets{g}, nSubjAll(g));
    fprintf('           (subject-level permutation requires >= 2 subjects).\n\n');
end
sheets         = sheets(~tooFew);
groupData      = groupData(~tooFew);
groupSubjNames = groupSubjNames(~tooFew);
nGroups        = numel(sheets);

if nGroups < 2
    error('Fewer than 2 groups with >= 2 subjects. Cannot run pairwise analysis.');
end

%% ── 3b. GUI: SELECT GROUPS FOR ANALYSIS ─────────────────────
rowH  = 30;
figW  = 420;
figH  = max(200, nGroups * rowH + 80);

selFig = figure('Name','Select Groups for Analysis', ...
                'NumberTitle','off', 'MenuBar','none', 'ToolBar','none', ...
                'Units','pixels', 'Position',[400 300 figW figH], ...
                'Resize','off', 'WindowStyle','modal');

uicontrol(selFig, 'Style','text', ...
          'String','Select groups to include in analysis:', ...
          'FontSize',10, 'FontWeight','bold', ...
          'HorizontalAlignment','left', ...
          'Units','pixels', 'Position',[15 figH-35 figW-20 22]);

chkHandles = gobjects(nGroups, 1);
for g = 1:nGroups
    yPos  = figH - 35 - g*rowH;
    nSub  = numel(groupData{g});
    nMin  = min(cellfun(@numel, groupData{g}));
    label = sprintf('%s   [%d subj, min %d samples/subj]', sheets{g}, nSub, nMin);
    chkHandles(g) = uicontrol(selFig, 'Style','checkbox', ...
                              'String', label, 'Value', 1, ...
                              'FontSize', 9, 'Units','pixels', ...
                              'Position', [15 yPos figW-30 24]);
end

uicontrol(selFig, 'Style','pushbutton', 'String','OK  -  Run Analysis', ...
          'FontSize',10, 'FontWeight','bold', ...
          'Units','pixels', 'Position',[figW/2-90 10 180 30], ...
          'Callback', @(~,~) uiresume(selFig));

uiwait(selFig);

if ~ishandle(selFig)
    error('Group selection window was closed. Script aborted.');
end
selectedMask = arrayfun(@(h) logical(h.Value), chkHandles);
close(selFig);

if sum(selectedMask) < 2
    error('At least 2 groups must be selected for pairwise K-S analysis.');
end

sheets         = sheets(selectedMask);
groupData      = groupData(selectedMask);
groupSubjNames = groupSubjNames(selectedMask);
nGroups        = numel(sheets);

fprintf('Groups selected for analysis (%d): %s\n\n', nGroups, strjoin(sheets, ', '));

%% ── 4. SAMPLE COUNTS + SMALLEST SUBJECT ─────────────────────
fprintf('%-30s  %-25s  %10s\n', 'Group', 'Subject', 'N_Samples');
fprintf('%s\n', repmat('-', 1, 70));

groupNSubjs = zeros(nGroups, 1);
groupMinN   = zeros(nGroups, 1);
for g = 1:nGroups
    counts = cellfun(@numel, groupData{g});
    for s = 1:numel(counts)
        fprintf('  %-28s  %-25s  %10d\n', sheets{g}, groupSubjNames{g}{s}, counts(s));
    end
    groupNSubjs(g) = numel(counts);
    groupMinN(g)   = min(counts);
    fprintf('  >>> Group "%s": %d subject(s), smallest subject N = %d\n\n', ...
            sheets{g}, groupNSubjs(g), groupMinN(g));
end
minSubjN = min(groupMinN);
fprintf('  Smallest per-subject sample count across selected groups: %d\n\n', minSubjN);

%% ── 5. USER DEFINES N_SAMP, N_ITER, N_PERM ──────────────────
dlgPrompts = { ...
    sprintf(['Samples drawn from EACH subject per iteration (N_SAMP).\n' ...
             'Drawn without replacement; identical for every subject in every group.\n' ...
             'Allowed range: 1 to %d (smallest per-subject sample count).'], minSubjN), ...
    'Number of subsampling iterations on the observed data (N_ITER):', ...
    sprintf(['Number of permutations (N_PERM).\n' ...
             'Subjects are shuffled between the two groups of each pair, and\n' ...
             'N_SAMP observations are drawn per subject exactly as above.'])};
dlgDefaults = {num2str(minSubjN), '10000', '10000'};

while true
    answer = inputdlg(dlgPrompts, 'Analysis Parameters', [1 70], dlgDefaults);
    if isempty(answer), error('User cancelled. Script aborted.'); end
    vals   = str2double(answer);
    N_SAMP = vals(1);
    N_ITER = vals(2);
    N_PERM = vals(3);

    msgs = {};
    if ~is_pos_int(N_SAMP) || N_SAMP > minSubjN
        msgs{end+1} = sprintf('N_SAMP must be a whole number from 1 to %d.', minSubjN); %#ok<SAGROW>
    end
    if ~is_pos_int(N_ITER), msgs{end+1} = 'N_ITER must be a positive whole number.'; end %#ok<SAGROW>
    if ~is_pos_int(N_PERM), msgs{end+1} = 'N_PERM must be a positive whole number.'; end %#ok<SAGROW>
    if isempty(msgs), break; end

    uiwait(errordlg(strjoin(msgs, newline), 'Invalid Parameters', 'modal'));
    dlgDefaults = answer(:)';
end

fprintf('Parameters confirmed:\n');
fprintf('  N_SAMP = %d  (per subject per iteration, without replacement)\n', N_SAMP);
fprintf('  N_ITER = %d\n  N_PERM = %d\n', N_ITER, N_PERM);
nFull = sum(cellfun(@(G) sum(cellfun(@numel, G) == N_SAMP), groupData));
if nFull > 0
    fprintf(['  Note: %d subject(s) have exactly N_SAMP observations and contribute\n' ...
             '        all of their data in every iteration.\n'], nFull);
end
fprintf('\n');

%% ── 5b. USER DEFINES CDF GRID ───────────────────────────────
allSubjVecs = vertcat(groupData{:});
allVals     = vertcat(allSubjVecs{:});
dataMin     = min(allVals);
dataMax     = max(allVals);

rawW = (dataMax - dataMin) / 100;
if rawW <= 0, rawW = 1; end
mag   = 10^floor(log10(rawW));
steps = [1 2 2.5 5 10] * mag;
defW  = steps(find(steps >= rawW, 1));
defStart = floor(dataMin / defW) * defW;
defEnd   = ceil(dataMax / defW) * defW;

gridPrompts = { ...
    sprintf(['Mean CDF grid: START value.\n' ...
             'Observed data range: %g to %g.'], dataMin, dataMax), ...
    'Mean CDF grid: END value (values above this are not shown on the CDF):', ...
    'Mean CDF grid: bin WIDTH:'};
gridDefaults = {num2str(defStart), num2str(defEnd), num2str(defW)};

while true
    answer = inputdlg(gridPrompts, 'CDF Grid', [1 60], gridDefaults);
    if isempty(answer), error('User cancelled. Script aborted.'); end
    gv = str2double(answer);
    BIN_START = gv(1); BIN_END = gv(2); BIN_WIDTH = gv(3);
    ok = all(isfinite(gv)) && BIN_WIDTH > 0 && BIN_END > BIN_START;
    if ok
        nSteps = round((BIN_END - BIN_START) / BIN_WIDTH);
        ok = (nSteps + 1) <= MAX_BINS;
    end
    if ok, break; end
    uiwait(errordlg(sprintf(['Grid must have END > START, WIDTH > 0, and at most ' ...
                             '%d points.'], MAX_BINS), 'Invalid CDF Grid', 'modal'));
    gridDefaults = answer(:)';
end
binEdges = BIN_START + (0:nSteps)' * BIN_WIDTH;
nBins    = numel(binEdges);
fprintf('CDF grid: %g to %g, width %g (%d points)\n\n', ...
        binEdges(1), binEdges(end), BIN_WIDTH, nBins);

%% ── 6. PAIRWISE LIST ─────────────────────────────────────────
pairs     = nchoosek(1:nGroups, 2);
nPairs    = size(pairs, 1);
pairNames = cell(nPairs, 1);
nDistinct = zeros(nPairs, 1);   % number of distinct subject reassignments
for p = 1:nPairs
    g1 = pairs(p,1);  g2 = pairs(p,2);
    pairNames{p} = sprintf('%s vs %s', sheets{g1}, sheets{g2});
    n1 = groupNSubjs(g1);  n2 = groupNSubjs(g2);
    nDistinct(p) = round(exp(gammaln(n1+n2+1) - gammaln(n1+1) - gammaln(n2+1)));
end
fprintf('%d pairwise comparison(s).\n', nPairs);
for p = find(nDistinct < 1000)'
    fprintf(['  Note: "%s" has only %d distinct ways to reassign subjects;\n' ...
             '        the permutation null is correspondingly coarse.\n'], ...
             pairNames{p}, nDistinct(p));
end
fprintf('\n');

pairColors = zeros(nPairs, 3);
for p = 1:nPairs
    pairColors(p, :) = hsv2rgb([(p - 1) / nPairs, 0.85, 0.88]);
end

%% ── 7. BALANCED SUBSAMPLING LOOP (OBSERVED DATA) ────────────
obsD    = zeros(N_ITER, nPairs);
cdfIter = cell(nGroups, 1);
for g = 1:nGroups
    cdfIter{g} = zeros(nBins, N_ITER);
end

fprintf('Running %d balanced subsampling iterations ...\n', N_ITER);
wb  = waitbar(0, 'Balanced subsampling (observed data)...', 'Name', 'K-S Analysis');
upd = max(1, round(N_ITER / 200));

for iter = 1:N_ITER
    pooled = cell(nGroups, 1);
    for g = 1:nGroups
        pooled{g} = draw_balanced(groupData{g}, N_SAMP);
        cdfIter{g}(:, iter) = binned_cdf(pooled{g}, binEdges);
    end
    for p = 1:nPairs
        obsD(iter, p) = ks2_stat(pooled{pairs(p,1)}, pooled{pairs(p,2)});
    end
    if mod(iter, upd) == 0 && isgraphics(wb)
        waitbar(iter / N_ITER, wb, sprintf('Iteration %d / %d', iter, N_ITER));
    end
end
if isgraphics(wb), close(wb); end
fprintf('  Done.\n\n');

%% ── 8. SUBJECT-LEVEL PERMUTATION LOOP (NULL) ────────────────
nullD = zeros(N_PERM, nPairs);

fprintf('Running %d permutations per pair ...\n', N_PERM);
wb    = waitbar(0, 'Permuting subjects...', 'Name', 'K-S Permutation');
total = nPairs * N_PERM;
upd   = max(1, round(total / 200));

for p = 1:nPairs
    g1 = pairs(p,1);  g2 = pairs(p,2);
    allSubs = [groupData{g1}; groupData{g2}];
    n1      = numel(groupData{g1});
    nTot    = numel(allSubs);

    for perm = 1:N_PERM
        order = randperm(nTot);
        pool1 = draw_balanced(allSubs(order(1:n1)),     N_SAMP);
        pool2 = draw_balanced(allSubs(order(n1+1:end)), N_SAMP);
        nullD(perm, p) = ks2_stat(pool1, pool2);

        step = (p - 1) * N_PERM + perm;
        if mod(step, upd) == 0 && isgraphics(wb)
            waitbar(step / total, wb, sprintf('Pair %d / %d,  permutation %d / %d', ...
                    p, nPairs, perm, N_PERM));
        end
    end
end
if isgraphics(wb), close(wb); end
fprintf('  Done.\n\n');

%% ── 8b. GLOBAL PERMUTATION FOR WESTFALL-YOUNG CORRECTION ────
% Subjects from ALL selected groups are shuffled together into
% pseudo-groups of the original sizes, so every pair's D comes from the
% same permutation. With only 2 groups this equals the pairwise null.
if nPairs > 1
    globalD  = zeros(N_PERM, nPairs);
    allSubsG = vertcat(groupData{:});
    nTotG    = numel(allSubsG);
    grpEnd   = cumsum(groupNSubjs(:))';
    grpStart = [1, grpEnd(1:end-1) + 1];

    fprintf('Running %d global permutations (Westfall-Young correction) ...\n', N_PERM);
    wb  = waitbar(0, 'Global permutation (Westfall-Young)...', 'Name', 'K-S Permutation');
    upd = max(1, round(N_PERM / 200));
    for perm = 1:N_PERM
        order  = randperm(nTotG);
        pooled = cell(nGroups, 1);
        for g = 1:nGroups
            pooled{g} = draw_balanced(allSubsG(order(grpStart(g):grpEnd(g))), N_SAMP);
        end
        for p = 1:nPairs
            globalD(perm, p) = ks2_stat(pooled{pairs(p,1)}, pooled{pairs(p,2)});
        end
        if mod(perm, upd) == 0 && isgraphics(wb)
            waitbar(perm / N_PERM, wb, sprintf('Global permutation %d / %d', perm, N_PERM));
        end
    end
    if isgraphics(wb), close(wb); end
    fprintf('  Done.\n\n');
end

%% ── 9. SUMMARY STATISTICS ───────────────────────────────────
qLevels = [ALPHA/2, 1 - ALPHA/2];

obsMeanD = mean(obsD, 1);
obsMedD  = median(obsD, 1);
obsSDD   = std(obsD, 0, 1);
obsCI    = quantile(obsD, qLevels, 1);        % 2 x nPairs

nullMeanD = mean(nullD, 1);
nullMedD  = median(nullD, 1);
nullCrit  = quantile(nullD, 1 - ALPHA, 1);    % one-sided critical D

% Empirical permutation p-value: how often the null D is at least as
% large as the mean observed D (+1 correction so p is never exactly 0).
pPerm = (1 + sum(nullD >= obsMeanD, 1)) ./ (1 + N_PERM);
[pBH,   sigBH]   = bh_correct(pPerm(:), ALPHA);
[pHolm, sigHolm] = holm_correct(pPerm(:), ALPHA);
if nPairs > 1
    pWY = westfall_young_minp(globalD, obsMeanD);
else
    pWY = pPerm(:);   % a single comparison needs no correction
end
sigWY = pWY <= ALPHA;
switch upper(PLOT_CORRECTION)
    case 'HOLM', sigPlot = sigHolm; plotLabel = 'Holm-corrected';
    case 'WY',   sigPlot = sigWY;   plotLabel = 'Westfall-Young-corrected';
    otherwise,   sigPlot = sigBH;   plotLabel = 'BH-corrected';
end

% Mean CDF per group across observed iterations
cdfStats = cell(nGroups, 1);
for g = 1:nGroups
    M = cdfIter{g};
    cdfStats{g}.mean = mean(M, 2);
    cdfStats{g}.sd   = std(M, 0, 2);
    cdfStats{g}.ci   = quantile(M, qLevels, 2);   % nBins x 2
end

% Console summary
ciPct = (1 - ALPHA) * 100;
fprintf('K-S summary (observed = balanced subsamples, null = subject permutation)\n');
fprintf('%-40s %10s %22s %10s %12s %14s %10s %10s %10s\n', 'Comparison', 'Obs mean D', ...
        sprintf('Obs %.0f%% interval', ciPct), 'Null mean', ...
        sprintf('Null %.0fth pct', (1-ALPHA)*100), 'p_uncorrected', 'p_BH', 'p_Holm', 'p_WY');
fprintf('%s\n', repmat('-', 1, 152));
for p = 1:nPairs
    fprintf('%-40s %10.4f   [%7.4f - %7.4f] %10.4f %12.4f %13.4g%-1s %9.4g%-1s %9.4g%-1s %9.4g%-1s\n', ...
            pairNames{p}, obsMeanD(p), obsCI(1,p), obsCI(2,p), nullMeanD(p), ...
            nullCrit(p), pPerm(p), ifelse_str(pPerm(p) <= ALPHA, '*', ''), ...
            pBH(p),   ifelse_str(sigBH(p),   '*', ''), ...
            pHolm(p), ifelse_str(sigHolm(p), '*', ''), ...
            pWY(p),   ifelse_str(sigWY(p),   '*', ''));
end
fprintf('BH = Benjamini-Hochberg (FDR); Holm = Holm-Bonferroni (FWER);\n');
fprintf('WY = Westfall-Young max-D / min-P permutation correction (FWER).\n');
fprintf('* = p <= %.2f\n', ALPHA);
fprintf('\n');

%% ── 10. EXPORT TO EXCEL ─────────────────────────────────────
fprintf('Exporting results to %s ...\n', OUTPUT_FILE);

% Start from a fresh workbook so no sheets from earlier runs remain
if isfile(OUTPUT_FILE)
    delete(OUTPUT_FILE);
    if isfile(OUTPUT_FILE)
        error('Cannot overwrite %s. Close it in Excel and re-run.', OUTPUT_FILE);
    end
end

% ---- Parameters ---------------------------------------------
paramSheet = { ...
    'Parameter',                                   'Value';
    'Input file',                                  EXCEL_FILE;
    'Groups analysed',                             strjoin(sheets, ', ');
    'N_SAMP (observations per subject per draw)',  N_SAMP;
    'Within-subject sampling',                     'Without replacement';
    'Smallest per-subject sample count',           minSubjN;
    'N_ITER (observed subsampling iterations)',    N_ITER;
    'N_PERM (permutations per pair)',              N_PERM;
    'Permutation scheme',                          'Subjects reassigned between the two groups of each pair';
    'p-value',                                     '(1 + #{null D >= mean observed D}) / (1 + N_PERM)';
    'Correction 1: Benjamini-Hochberg',            'False discovery rate, across pairs';
    'Correction 2: Holm-Bonferroni',               'Family-wise error rate, across pairs';
    'Correction 3: Westfall-Young',                'Family-wise error rate; single-step min-P over D, subjects permuted across all groups';
    'Global permutations (Westfall-Young)',        ifelse_num(nPairs > 1, N_PERM, 0);
    'Correction marked on summary plot',           plotLabel;
    'ALPHA',                                       ALPHA;
    'Random seed (0 = unseeded)',                  SEED;
    'CDF grid start',                              binEdges(1);
    'CDF grid end',                                binEdges(end);
    'CDF grid width',                              BIN_WIDTH;
    'Run date',                                    char(datetime('now'))};
writecell(paramSheet, OUTPUT_FILE, 'Sheet', 'Parameters');

% ---- Subject statistics -------------------------------------
statHeaders = {'Group','Subject','N_Samples','Fraction_Drawn_Per_Iter', ...
               'Min','Max','Mean','Median','Std_Dev','SEM','CV_pct'};
statRows = cell(0, numel(statHeaders));
for g = 1:nGroups
    for s = 1:numel(groupData{g})
        x  = groupData{g}{s};
        n  = numel(x);
        mn = mean(x);
        sd = std(x);
        statRows(end+1, :) = {sheets{g}, groupSubjNames{g}{s}, n, N_SAMP / n, ...
                              min(x), max(x), mn, median(x), sd, sd / sqrt(n), ...
                              100 * sd / mn}; %#ok<SAGROW>
    end
end
writecell([statHeaders; statRows], OUTPUT_FILE, 'Sheet', 'Subject_Statistics');

% ---- K-S summary --------------------------------------------
sumHeaders = {'Comparison','Group_1','Group_2','N_Subjects_G1','N_Subjects_G2', ...
              'N_Distinct_Subject_Assignments', ...
              'Obs_Mean_D','Obs_Median_D','Obs_SD_D', ...
              sprintf('Obs_D_Lo_%.0fPct', ciPct), sprintf('Obs_D_Hi_%.0fPct', ciPct), ...
              'Null_Mean_D','Null_Median_D', ...
              sprintf('Null_Critical_D_%.0fthPct', (1-ALPHA)*100), ...
              'p_uncorrected','Sig_uncorrected','p_BH_corrected','Sig_BH_corrected', ...
              'p_Holm_corrected','Sig_Holm_corrected', ...
              'p_WestfallYoung_corrected','Sig_WestfallYoung_corrected'};
sumRows = cell(nPairs, numel(sumHeaders));
for p = 1:nPairs
    g1 = pairs(p,1);  g2 = pairs(p,2);
    sumRows(p, :) = {pairNames{p}, sheets{g1}, sheets{g2}, ...
                     groupNSubjs(g1), groupNSubjs(g2), nDistinct(p), ...
                     obsMeanD(p), obsMedD(p), obsSDD(p), obsCI(1,p), obsCI(2,p), ...
                     nullMeanD(p), nullMedD(p), nullCrit(p), ...
                     pPerm(p), ifelse_str(pPerm(p) <= ALPHA, '*', ''), ...
                     pBH(p),   ifelse_str(sigBH(p),   '*', ''), ...
                     pHolm(p), ifelse_str(sigHolm(p), '*', ''), ...
                     pWY(p),   ifelse_str(sigWY(p),   '*', '')};
end
writecell([sumHeaders; sumRows], OUTPUT_FILE, 'Sheet', 'KS_Summary');

% ---- Full observed D matrix [N_ITER x nPairs] ----------------
writecell([[{'Iteration'}, pairNames']; num2cell([(1:N_ITER)', obsD])], ...
          OUTPUT_FILE, 'Sheet', 'Observed_D_Iterations');

% ---- Full permutation null D matrix [N_PERM x nPairs] --------
writecell([[{'Permutation'}, pairNames']; num2cell([(1:N_PERM)', nullD])], ...
          OUTPUT_FILE, 'Sheet', 'Permutation_D_Null');

% ---- Global permutation D matrix (Westfall-Young) -------------
if nPairs > 1
    writecell([[{'Global_Permutation'}, pairNames']; num2cell([(1:N_PERM)', globalD])], ...
              OUTPUT_FILE, 'Sheet', 'Global_Permutation_D_Null');
end

% ---- Mean CDF for each pairwise comparison (observed data) ---
% One block per pair, side by side, separated by a blank column.
blockW  = 10;
sheetW  = nPairs * (blockW + 1) - 1;
cdfSheet = cell(nBins + 2, sheetW);
for p = 1:nPairs
    g1 = pairs(p,1);  g2 = pairs(p,2);
    n1 = sheets{g1};  n2 = sheets{g2};
    c0 = (p - 1) * (blockW + 1);
    cols = c0 + (1:blockW);

    cdfSheet{1, c0 + 1} = pairNames{p};
    cdfSheet(2, cols) = {'Bin_End', ...
        [n1 '_MeanCDF'], [n1 '_SD'], ...
        sprintf('%s_Lo_%.0fPct', n1, ciPct), sprintf('%s_Hi_%.0fPct', n1, ciPct), ...
        [n2 '_MeanCDF'], [n2 '_SD'], ...
        sprintf('%s_Lo_%.0fPct', n2, ciPct), sprintf('%s_Hi_%.0fPct', n2, ciPct), ...
        'MeanCDF_Diff_G1_minus_G2'};
    block = [binEdges, ...
             cdfStats{g1}.mean, cdfStats{g1}.sd, cdfStats{g1}.ci, ...
             cdfStats{g2}.mean, cdfStats{g2}.sd, cdfStats{g2}.ci, ...
             cdfStats{g1}.mean - cdfStats{g2}.mean];
    cdfSheet(3:end, cols) = num2cell(block);
end
writecell(cdfSheet, OUTPUT_FILE, 'Sheet', 'Mean_CDF_by_Pair');

fprintf('  Export complete.\n\n');

%% ── 11. FIGURES ─────────────────────────────────────────────
colG1 = [0.000 0.447 0.741];
colG2 = [0.850 0.325 0.098];

% ---- 11a. Observed D vs permutation null, one figure per pair --
figPair = gobjects(nPairs, 1);
for p = 1:nPairs
    figPair(p) = figure('Name', sprintf('Observed vs Null D - %s', pairNames{p}), ...
                        'Color','w', 'Units','normalized', 'Position',[0.1 0.15 0.7 0.6]);
    ax = axes(figPair(p));  hold(ax, 'on');

    hiD = max([obsD(:,p); nullD(:,p)]);
    if hiD <= 0, hiD = 1; end
    edges = linspace(0, hiD * 1.05, 61);

    h1 = histogram(ax, nullD(:,p), edges, 'Normalization','probability', ...
                   'FaceColor',[0.6 0.6 0.6], 'EdgeColor','none', 'FaceAlpha',0.6);
    h2 = histogram(ax, obsD(:,p), edges, 'Normalization','probability', ...
                   'FaceColor',pairColors(p,:), 'EdgeColor','none', 'FaceAlpha',0.6);
    h3 = xline(ax, obsMeanD(p), '-r', 'LineWidth', 2);
    h4 = xline(ax, nullCrit(p), '--k', 'LineWidth', 1.5);

    xlabel(ax, 'K-S statistic D');
    ylabel(ax, 'Proportion of iterations');
    legend([h1 h2 h3 h4], { ...
        sprintf('Permutation null (N_{perm} = %d)', N_PERM), ...
        sprintf('Observed subsamples (N_{iter} = %d)', N_ITER), ...
        sprintf('Observed mean D = %.4f', obsMeanD(p)), ...
        sprintf('Null %.0fth percentile = %.4f', (1-ALPHA)*100, nullCrit(p))}, ...
        'Location','best');
    title(ax, {pairNames{p}, ...
               sprintf('p uncorrected = %.4g,   BH = %.4g,   Holm = %.4g,   Westfall-Young = %.4g', ...
               pPerm(p), pBH(p), pHolm(p), pWY(p))}, 'Interpreter','none');
    box(ax, 'off');  grid(ax, 'on');
    drawnow;
end

% ---- 11b. Mean CDF curves, one figure per pair -----------------
figCDF = gobjects(nPairs, 1);
for p = 1:nPairs
    g1 = pairs(p,1);  g2 = pairs(p,2);
    figCDF(p) = figure('Name', sprintf('Mean CDF - %s', pairNames{p}), ...
                       'Color','w', 'Units','normalized', 'Position',[0.1 0.1 0.7 0.65]);
    ax = axes(figCDF(p));  hold(ax, 'on');
    be = binEdges;

    f1 = fill(ax, [be; flipud(be)], [cdfStats{g1}.ci(:,1); flipud(cdfStats{g1}.ci(:,2))], ...
              colG1, 'FaceAlpha',0.18, 'EdgeColor','none');
    l1 = plot(ax, be, cdfStats{g1}.mean, '-',  'Color', colG1, 'LineWidth', 2.2);
    f2 = fill(ax, [be; flipud(be)], [cdfStats{g2}.ci(:,1); flipud(cdfStats{g2}.ci(:,2))], ...
              colG2, 'FaceAlpha',0.18, 'EdgeColor','none');
    l2 = plot(ax, be, cdfStats{g2}.mean, '--', 'Color', colG2, 'LineWidth', 2.2);

    xlabel(ax, sprintf('Value (grid width = %g)', BIN_WIDTH));
    ylabel(ax, 'Cumulative probability');
    ylim(ax, [0 1]);
    legend([l1 f1 l2 f2], { ...
        sprintf('%s mean CDF', sheets{g1}), sprintf('%s %.0f%% interval', sheets{g1}, ciPct), ...
        sprintf('%s mean CDF', sheets{g2}), sprintf('%s %.0f%% interval', sheets{g2}, ciPct)}, ...
        'Interpreter','none', 'Location','best');
    title(ax, sprintf('Mean cumulative distribution (N_SAMP = %d per subject, %d iterations) - %s', ...
          N_SAMP, N_ITER, pairNames{p}), 'Interpreter','none');
    box(ax, 'off');  grid(ax, 'on');
    drawnow;
end

% ---- 11c. Summary: mean observed D vs null critical D ----------
figSum = figure('Name','Mean observed D vs null', 'Color','w', ...
                'Units','normalized', 'Position',[0.1 0.2 0.75 0.55]);
ax = axes(figSum);  hold(ax, 'on');
for p = 1:nPairs
    bar(ax, p, obsMeanD(p), 'FaceColor', pairColors(p,:), 'FaceAlpha', 0.85);
end
hErr  = errorbar(ax, 1:nPairs, obsMeanD, obsMeanD - obsCI(1,:), obsCI(2,:) - obsMeanD, ...
                 'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',8);
hCrit = plot(ax, 1:nPairs, nullCrit, 'kd', 'MarkerFaceColor','w', 'MarkerSize',9, 'LineWidth',1.5);
for p = find(sigPlot(:))'
    text(ax, p, obsCI(2,p), '*', 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontSize',18, 'FontWeight','bold');
end
set(ax, 'XTick', 1:nPairs, 'XTickLabel', pairNames, 'XTickLabelRotation', 30, ...
        'TickLabelInterpreter','none');
ylabel(ax, 'K-S statistic D');
legend([hErr hCrit], {sprintf('Observed %.0f%% interval', ciPct), ...
        sprintf('Null %.0fth percentile', (1-ALPHA)*100)}, 'Location','best');
title(ax, sprintf('Mean observed D per comparison  (* = %s p <= %.2f)', plotLabel, ALPHA));
box(ax, 'off');  grid(ax, 'on');
drawnow;

%% ── 12. SAVE FIGURES ────────────────────────────────────────
fprintf('Saving figures ...\n');
for p = 1:nPairs
    tag = regexprep(pairNames{p}, '[^\w]+', '_');
    save_fig(figPair(p), fullfile(fpath, [baseName OUT_TAG '_ObsVsNull_' tag '.png']));
    save_fig(figCDF(p),  fullfile(fpath, [baseName OUT_TAG '_MeanCDF_'   tag '.png']));
end
save_fig(figSum, fullfile(fpath, [baseName OUT_TAG '_Summary.png']));

close all;
fprintf('\nDone.\n');


%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================

function pooled = draw_balanced(subjList, k)
% Draw k observations WITHOUT replacement from every subject in
% subjList (cell array of vectors) and stack them into one column.
    nS     = numel(subjList);
    pooled = zeros(nS * k, 1);
    for s = 1:nS
        x = subjList{s};
        pooled((s-1)*k + 1 : s*k) = x(randperm(numel(x), k));
    end
end

function D = ks2_stat(x, y)
% Two-sample K-S statistic D = max |F_x(t) - F_y(t)|, evaluated at
% every distinct pooled value. Equivalent to the statistic returned
% by kstest2 (two-sided), without computing the asymptotic p-value.
    x = x(:);  y = y(:);
    u  = unique([x; y]);
    c1 = histcounts(x, [u; Inf]);   % counts of x at each distinct value
    c2 = histcounts(y, [u; Inf]);
    D  = max(abs(cumsum(c1) / numel(x) - cumsum(c2) / numel(y)));
end

function F = binned_cdf(x, edges)
% Empirical CDF P(X <= edge) at each grid edge.
    b      = discretize(x, [-Inf; edges; Inf], 'IncludedEdge', 'right');
    counts = accumarray(b, 1, [numel(edges) + 1, 1]);
    F      = cumsum(counts(1:numel(edges))) / numel(x);
end

function [p_adj, sig_mask] = bh_correct(pvals, q)
% Benjamini-Hochberg FDR-adjusted p-values (monotone, capped at 1).
    m        = numel(pvals);
    [ps, si] = sort(pvals(:), 'ascend');
    p_mono   = ps .* m ./ (1:m)';
    for i = m-1 : -1 : 1
        p_mono(i) = min(p_mono(i), p_mono(i+1));
    end
    p_mono    = min(p_mono, 1);
    p_adj     = zeros(m, 1);
    p_adj(si) = p_mono;
    sig_mask  = p_adj <= q;
end

function [p_adj, sig_mask] = holm_correct(pvals, a)
% Holm-Bonferroni step-down adjusted p-values (monotone, capped at 1).
    m        = numel(pvals);
    [ps, si] = sort(pvals(:), 'ascend');
    p_step   = ps .* (m - (1:m)' + 1);
    for i = 2:m
        p_step(i) = max(p_step(i), p_step(i-1));
    end
    p_step    = min(p_step, 1);
    p_adj     = zeros(m, 1);
    p_adj(si) = p_step;
    sig_mask  = p_adj <= a;
end

function p_adj = westfall_young_minp(nullMat, obsStat)
% Westfall-Young single-step min-P adjusted p-values.
%   nullMat : [nPerm x nPairs] statistics from permutations that shuffle
%             subjects across ALL groups at once.
%   obsStat : [1 x nPairs] observed statistic (here, mean observed D).
% Every value is converted to a null p-value within its own pair, using
% the same (1 + count) / (1 + nPerm) convention as the observed p, so
% pairs with differently shaped null distributions are on one scale.
% The adjusted p is the fraction of permutations whose smallest p
% across pairs is <= the pair's observed p.
    [nPerm, nP] = size(nullMat);
    pNull = zeros(nPerm, nP);
    pObs  = zeros(nP, 1);
    for k = 1:nP
        d = nullMat(:, k);
        [~, ~, ic] = unique(d);
        cnt  = accumarray(ic, 1);
        geq  = flipud(cumsum(flipud(cnt)));   % #{null >= each distinct value}
        pNull(:, k) = geq(ic) / (nPerm + 1);  % itself plus others >= it
        pObs(k)     = (1 + sum(d >= obsStat(k))) / (nPerm + 1);
    end
    minP  = min(pNull, [], 2);
    p_adj = (1 + sum(minP <= pObs', 1))' / (nPerm + 1);
    p_adj = min(p_adj, 1);
end

function v = ifelse_num(condition, yes_val, no_val)
    if condition, v = yes_val; else, v = no_val; end
end

function tf = is_blank_cell(x)
% True for empty / missing / whitespace-only cells from readcell.
    tf = isempty(x) || isa(x, 'missing') || ...
         ((ischar(x) || isstring(x)) && strlength(strtrim(string(x))) == 0);
end

function tf = is_pos_int(v)
    tf = isscalar(v) && isfinite(v) && v >= 1 && v == round(v);
end

function s = ifelse_str(condition, yes_val, no_val)
    if condition, s = yes_val; else, s = no_val; end
end

function save_fig(fig, fileName)
    try
        exportgraphics(fig, fileName, 'Resolution', 300);
        fprintf('  Saved %s\n', fileName);
    catch err
        fprintf('  Could not save %s (%s)\n', fileName, err.message);
    end
end

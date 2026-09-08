function R = extract_stim_peristim_core(z, fs, onset, off, sumimg, Sn, frame_px, out_dir, run_id, params)
%EXTRACT_STIM_PERISTIM_CORE  Peri-stimulus analysis + figures (format-agnostic).
%
%   R = extract_stim_peristim_core(z, fs, onset, off, sumimg, Sn, frame_px, ...
%                                  out_dir, run_id, params)
%
% Shared engine used by the format-specific loaders
% (extract_stim_triggered_analysis = moco_frr, extract_stim_ext_analysis = _ext).
%
% INPUTS
%   z        : frames x neurons, z-scored activity (already normalized).
%   fs       : trace frame rate (Hz).
%   onset    : vector of poke-onset frame indices IN TRACE FRAME SPACE (sorted).
%   off      : struct with .pre_s / .post_s = mean pre/post offsets (s) for
%              markers ([] or NaN to omit); .ef_lbl label for the t=0 field.
%   sumimg   : summary image for FOV context ([] to omit).
%   Sn       : normalized sparse spatial matrix (n_pix x n_neurons) for the FOV
%              response map, or [] to skip the spatial figure.
%   frame_px : square FOV size (used only if Sn provided).
%   out_dir  : output folder (created if missing).
%   run_id   : label for titles/filenames.
%   params   : struct (peri_win, base_win, resp_win, z_thresh, resp_test,
%              alpha, max_show). See loaders for defaults.

peri_win  = getp(params, 'peri_win', [-2 3]);
base_win  = getp(params, 'base_win', [-2 0]);
resp_win  = getp(params, 'resp_win', [0 2]);
z_thresh  = getp(params, 'z_thresh', 2);
resp_test = getp(params, 'resp_test', 'signrank');
alpha     = getp(params, 'alpha', 0.05);
max_show  = getp(params, 'max_show', 2000);
ef_lbl    = getp(off, 'ef_lbl', 'onset');
pre_off_s  = getp(off, 'pre_s', NaN);
post_off_s = getp(off, 'post_s', NaN);
col_exc = [0.85 0.20 0.15];
col_inh = [0.10 0.35 0.85];

if exist(out_dir, 'dir') ~= 7; mkdir(out_dir); end
[n_frames, n_neurons] = size(z);

% --- Peri-stimulus window (frames) --------------------------------------
pk = round(peri_win(1)*fs) : round(peri_win(2)*fs);
nP = numel(pk);
t_peri = pk(:) / fs;
idx_base = t_peri >= base_win(1) & t_peri < base_win(2);
idx_resp = t_peri >= resp_win(1) & t_peri < resp_win(2);
if ~any(idx_base) || ~any(idx_resp)
    error('extract_stim_peristim_core:EmptyWindow', ...
        'base_win/resp_win fall outside peri_win; check parameters.');
end

onset = sort(round(double(onset(:))));
n_pokes_total = numel(onset);
valid = (onset + pk(1) >= 1) & (onset + pk(end) <= n_frames);
onset_v = onset(valid);
nT = numel(onset_v);
if nT == 0
    error('extract_stim_peristim_core:NoValidTrials', ...
        'No pokes have a full peri-window inside the recording.');
end
fprintf('  %d/%d pokes usable (peri-window in range)\n', nT, n_pokes_total);

% --- Accumulate STA + per-trial windows ---------------------------------
sta_sum     = zeros(n_neurons, nP);
pop_trials  = zeros(nT, nP);
frac_trials = zeros(nT, nP);
trial_base  = zeros(n_neurons, nT);
trial_resp  = zeros(n_neurons, nT);
for k = 1:nT
    seg = double(z(onset_v(k) + pk, :));       % nP x n_neurons
    sta_sum = sta_sum + seg.';
    pop_trials(k, :)  = mean(seg, 2).';
    frac_trials(k, :) = mean(seg > z_thresh, 2).';
    trial_base(:, k)  = mean(seg(idx_base, :), 1).';
    trial_resp(:, k)  = mean(seg(idx_resp, :), 1).';
end
sta = sta_sum / nT;
clear sta_sum

% --- Per-neuron responsiveness (paired pre vs post) ---------------------
fprintf('  Responsiveness test (%s) over %d neurons...\n', resp_test, n_neurons);
amp = mean(trial_resp - trial_base, 2);
pval = ones(n_neurons, 1);
if strcmpi(resp_test, 'ttest')
    d = trial_resp - trial_base;
    md = mean(d, 2); sd = std(d, 0, 2);
    tstat = md ./ max(sd ./ sqrt(nT), 1e-12);
    pval = 2 * tcdf(-abs(tstat), max(nT-1, 1));
    pval(sd < 1e-12) = 1;
else
    for nrn = 1:n_neurons
        d = trial_resp(nrn, :) - trial_base(nrn, :);
        if all(abs(d) < 1e-12)
            pval(nrn) = 1;
        else
            pval(nrn) = signrank(trial_resp(nrn, :), trial_base(nrn, :));
        end
    end
end
is_exc = pval < alpha & amp > 0;
is_inh = pval < alpha & amp < 0;
n_exc = sum(is_exc); n_inh = sum(is_inh);

pop_base = mean(trial_base, 1).';
pop_resp = mean(trial_resp, 1).';
if all(abs(pop_resp - pop_base) < 1e-12)
    p_pop = 1;
else
    p_pop = signrank(pop_resp, pop_base);
end
mean_pop_resp = mean(pop_resp - pop_base);

if n_pokes_total >= 2; iti = diff(onset) / fs; else; iti = []; end

% --- Colormaps + SEM -----------------------------------------------------
[cmap_inf, cmap_rb, cmap_vir] = stim_colormaps();
sem = @(x) std(x, 0, 1) ./ sqrt(max(size(x, 1), 1));

% ======================= FIGURE 1: population PSTH ======================
fprintf('  Plotting figures...\n');
figure('Position', [30 30 1500 560], 'Color', 'w');
sgtitle(sprintf('%s  |  Population peri-poke response (t=0 at %s, n=%d pokes)', ...
    run_id, ef_lbl, nT), 'FontSize', 13, 'FontWeight', 'bold');

subplot(1, 2, 1);
mu = mean(pop_trials, 1); se = sem(pop_trials);
fill([t_peri; flipud(t_peri)], [mu(:)-se(:); flipud(mu(:)+se(:))], ...
    col_exc, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
plot(t_peri, mu, 'Color', col_exc, 'LineWidth', 2);
add_event_markers(pre_off_s, post_off_s);
xlabel('Time from poke (s)'); ylabel('Population mean z-score');
title(sprintf('Population activity (mean \\pm SEM)\nresp vs base p = %.2e', p_pop));
box off;

subplot(1, 2, 2);
mu = mean(frac_trials, 1); se = sem(frac_trials);
fill([t_peri; flipud(t_peri)], [mu(:)-se(:); flipud(mu(:)+se(:))], ...
    [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
plot(t_peri, mu, 'k', 'LineWidth', 2);
add_event_markers(pre_off_s, post_off_s);
xlabel('Time from poke (s)'); ylabel(sprintf('Fraction active (z>%g)', z_thresh));
title('Population recruitment'); box off;

exportgraphics(gcf, fullfile(out_dir, 'fig1_stim_population.png'), 'Resolution', 200);

% ======================= FIGURE 2: per-neuron STA ======================
figure('Position', [30 30 1500 620], 'Color', 'w');
sgtitle(sprintf('%s  |  Stimulus-triggered average per neuron', run_id), ...
    'FontSize', 13, 'FontWeight', 'bold');

[~, ord] = sort(amp, 'descend');
if n_neurons > max_show
    ord_show = ord(round(linspace(1, n_neurons, max_show)));
else
    ord_show = ord;
end
subplot(1, 2, 1);
imagesc(peri_win, [1 numel(ord_show)], sta(ord_show, :));
caxis([-1 3]); colormap(gca, cmap_inf);
cb = colorbar; cb.Label.String = 'z-score';
hold on; xline(0, 'w--', 'LineWidth', 1);
xlabel('Time from poke (s)'); ylabel('Neuron (sorted by response amp.)');
title(sprintf('STA heatmap (%d neurons)', numel(ord_show)));

subplot(1, 2, 2);
leg = {};
if n_exc > 0; plot(t_peri, mean(sta(is_exc, :), 1), 'Color', col_exc, 'LineWidth', 2); hold on; leg{end+1} = sprintf('Excited (n=%d)', n_exc); end
if n_inh > 0; plot(t_peri, mean(sta(is_inh, :), 1), 'Color', col_inh, 'LineWidth', 2); hold on; leg{end+1} = sprintf('Inhibited (n=%d)', n_inh); end
plot(t_peri, mean(sta, 1), 'k', 'LineWidth', 1.5); leg{end+1} = 'All neurons';
add_event_markers(pre_off_s, post_off_s);
xlabel('Time from poke (s)'); ylabel('Mean STA (z-score)');
legend(leg, 'Location', 'best'); box off; title('Group-mean STA');

exportgraphics(gcf, fullfile(out_dir, 'fig2_neuron_sta_heatmap.png'), 'Resolution', 200);

% ======================= FIGURE 3: responsive cells ====================
figure('Position', [30 30 1500 500], 'Color', 'w');
sgtitle(sprintf('%s  |  Poke-responsive neurons (%s, alpha=%.3g)', ...
    run_id, resp_test, alpha), 'FontSize', 13, 'FontWeight', 'bold');

subplot(1, 3, 1);
n_ns = n_neurons - n_exc - n_inh;
b = bar([1 2 3], [n_exc n_inh n_ns], 0.6, 'EdgeColor', 'k'); hold on;
b.FaceColor = 'flat'; b.CData = [col_exc; col_inh; 0.7 0.7 0.7];
set(gca, 'XTick', 1:3, 'XTickLabel', {'Excited', 'Inhibited', 'n.s.'});
ylabel('# neurons');
text(1, double(n_exc), sprintf('%.1f%%', 100*n_exc/n_neurons), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(2, double(n_inh), sprintf('%.1f%%', 100*n_inh/n_neurons), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
title('Responsive cell counts'); box off;

subplot(1, 3, 2);
edges = linspace(prctile(amp, 1), prctile(amp, 99), 60);
histogram(amp, edges, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); hold on;
xline(0, 'k--', 'LineWidth', 1); xline(mean(amp), 'r-', 'LineWidth', 1.5);
xlabel('Response amplitude (post - pre, z)'); ylabel('# neurons');
title(sprintf('Amplitude distribution\nmean = %.4f', mean(amp))); box off;

subplot(1, 3, 3);
[fa, xa] = ecdf(pval);
plot(xa, fa, 'k', 'LineWidth', 2); hold on;
xline(alpha, 'r--', 'LineWidth', 1.5);
plot([0 1], [0 1], ':', 'Color', [0.6 0.6 0.6]);
xlabel('p-value'); ylabel('CDF'); title('p-value distribution'); box off; xlim([0 1]);

exportgraphics(gcf, fullfile(out_dir, 'fig3_responsive_cells.png'), 'Resolution', 200);

% ======================= FIGURE 4: spatial response ====================
if ~isempty(Sn)
    P4 = struct('Sn', Sn, 'amp', amp, 'isExc', is_exc, 'isInh', is_inh, ...
        'framePx', frame_px, 'sumimg', sumimg, ...
        'outFile', fullfile(out_dir, 'fig4_spatial_response.png'), 'runId', run_id);
    R_map_small = extract_stim_plot_spatial_fig(P4);
elseif ~isempty(sumimg)
    figure('Position', [30 30 700 620], 'Color', 'w');
    imagesc(sumimg); colormap(gca, 'gray'); axis image off;
    title(sprintf('%s  |  Summary image (no spatial footprints available)', run_id));
    exportgraphics(gcf, fullfile(out_dir, 'fig4_summary_image.png'), 'Resolution', 200);
    R_map_small = [];
else
    R_map_small = [];
end

% ======================= FIGURE 5: trial consistency ===================
figure('Position', [30 30 1400 520], 'Color', 'w');
sgtitle(sprintf('%s  |  Trial consistency & timing', run_id), 'FontSize', 13, 'FontWeight', 'bold');

subplot(1, 2, 1);
imagesc(peri_win, [1 nT], pop_trials);
colormap(gca, cmap_vir); cb = colorbar; cb.Label.String = 'Pop. mean z';
hold on; xline(0, 'w--', 'LineWidth', 1);
xlabel('Time from poke (s)'); ylabel('Poke # (trial)');
title(sprintf('Per-trial population activity (n=%d)', nT));

subplot(1, 2, 2);
if ~isempty(iti)
    histogram(iti, 'FaceColor', [0.4 0.4 0.7], 'EdgeColor', 'none');
    xlabel('Inter-poke interval (s)'); ylabel('Count');
    title(sprintf('Poke intervals (median = %.2f s, n = %d)', median(iti), numel(iti)));
else
    text(0.5, 0.5, 'Not enough pokes for ITI', 'HorizontalAlignment', 'center'); axis off;
end
box off;

exportgraphics(gcf, fullfile(out_dir, 'fig5_trial_consistency.png'), 'Resolution', 200);

% ======================= summary text ==================================
fid = fopen(fullfile(out_dir, 'stim_summary.txt'), 'w');
L = {};
L{end+1} = sprintf('%s | Stimulus-triggered (facial poke) analysis', run_id);
L{end+1} = repmat('=', 1, 78);
L{end+1} = sprintf('%-40s %s', 'Onset event field (t=0)', ef_lbl);
L{end+1} = sprintf('%-40s %d', 'Neurons', n_neurons);
L{end+1} = sprintf('%-40s %d', 'Frames', n_frames);
L{end+1} = sprintf('%-40s %.3f', 'Frame rate (Hz)', fs);
L{end+1} = sprintf('%-40s %d', 'Pokes detected', n_pokes_total);
L{end+1} = sprintf('%-40s %d', 'Pokes used (peri-window in range)', nT);
if ~isempty(iti)
    L{end+1} = sprintf('%-40s %.2f (%.2f-%.2f)', 'Inter-poke interval s (med, range)', median(iti), min(iti), max(iti));
end
L{end+1} = sprintf('%-40s [%.2f %.2f] s', 'Peri window', peri_win(1), peri_win(2));
L{end+1} = sprintf('%-40s [%.2f %.2f] s', 'Baseline window', base_win(1), base_win(2));
L{end+1} = sprintf('%-40s [%.2f %.2f] s', 'Response window', resp_win(1), resp_win(2));
L{end+1} = '';
L{end+1} = '--- Responsiveness ---';
L{end+1} = sprintf('%-40s %d (%.1f%%)', 'Excited neurons', n_exc, 100*n_exc/n_neurons);
L{end+1} = sprintf('%-40s %d (%.1f%%)', 'Inhibited neurons', n_inh, 100*n_inh/n_neurons);
L{end+1} = sprintf('%-40s %.4f', 'Mean response amp (post-pre, z)', mean(amp));
L{end+1} = sprintf('%-40s %.4f', 'Population resp (post-pre, z)', mean_pop_resp);
L{end+1} = sprintf('%-40s %.3e', 'Population resp vs base p', p_pop);
for k = 1:numel(L); fprintf('%s\n', L{k}); fprintf(fid, '%s\n', L{k}); end
fclose(fid);

% ======================= return struct =================================
R = struct();
R.run_id = run_id; R.event_field = ef_lbl; R.fs = fs;
R.n_neurons = n_neurons; R.n_pokes_total = n_pokes_total; R.n_trials = nT;
R.n_exc = n_exc; R.n_inh = n_inh;
R.frac_exc = n_exc / n_neurons; R.frac_inh = n_inh / n_neurons;
R.mean_amp = mean(amp); R.mean_pop_resp = mean_pop_resp; R.p_pop = p_pop;
if ~isempty(iti); R.iti_median = median(iti); else; R.iti_median = NaN; end
R.t_peri = t_peri; R.pop_psth = mean(pop_trials, 1).'; R.frac_psth = mean(frac_trials, 1).';
R.sta_all = mean(sta, 1).';
if n_exc > 0; R.sta_exc = mean(sta(is_exc, :), 1).'; else; R.sta_exc = nan(nP, 1); end
if n_inh > 0; R.sta_inh = mean(sta(is_inh, :), 1).'; else; R.sta_inh = nan(nP, 1); end
R.amp = amp; R.pval = pval; R.iti = iti;
if ~isempty(sumimg)
    R.sumimg_small = imresize(double(sumimg), [256 256]);
else
    R.sumimg_small = [];
end
R.resp_map_small = R_map_small;

end

% ======================= local helpers =================================
function v = getp(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function add_event_markers(pre_off_s, post_off_s)
xline(0, 'k-', 'LineWidth', 1.5);
if ~isempty(pre_off_s) && isfinite(pre_off_s); xline(pre_off_s, 'b:', 'LineWidth', 1); end
if ~isempty(post_off_s) && isfinite(post_off_s); xline(post_off_s, 'm:', 'LineWidth', 1); end
end

function [cmap_inf, cmap_rb, cmap_vir] = stim_colormaps()
t = linspace(0, 1, 256)';
inf_r = min(1, max(0, -0.0155 + t.*(5.27 + t.*(-14.35 + t.*(17.77 + t.*(-9.42 + t.*1.76))))));
inf_g = min(1, max(0, -0.0123 + t.*(-0.23 + t.*(7.77 + t.*(-17.05 + t.*(16.95 + t.*(-6.34)))))));
inf_b = min(1, max(0,  0.095 + t.*(4.00 + t.*(-15.24 + t.*(23.56 + t.*(-17.39 + t.*5.01))))));
cmap_inf = [inf_r inf_g inf_b];
half_n = 128;
cmap_rb = [[linspace(0, 1, half_n) ones(1, half_n)]' ...
           [linspace(0, 1, half_n) linspace(1, 0, half_n)]' ...
           [ones(1, half_n) linspace(1, 0, half_n)]'];
vr = min(1, max(0, 0.267 + t.*(-0.005 + t.*(3.07 + t.*(-6.36 + t.*(5.27 + t.*(-1.58)))))));
vg = min(1, max(0, 0.004 + t.*(1.42 + t.*(-0.15 + t.*(-1.95 + t.*(3.11 + t.*(-1.40)))))));
vb = min(1, max(0, 0.329 + t.*(1.44 + t.*(-5.10 + t.*(9.18 + t.*(-8.28 + t.*2.75))))));
cmap_vir = [vr vg vb];
end

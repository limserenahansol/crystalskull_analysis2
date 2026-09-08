function extract_stim_lab_meeting_figures(results, out_root, params)
%EXTRACT_STIM_LAB_MEETING_FIGURES  5 per-recording copies + 5 cohort averages.
%
%   extract_stim_lab_meeting_figures(results, out_root, params)
%
% Writes:
%   out_root/lab_meeting/per_recording/<label>/fig1..fig5
%   out_root/lab_meeting/cohort_average/fig1..fig5
%
% Fig4 is the poke-response map (spatial_weights x per-neuron amplitude).
% Fallback: EXTRACT summary image if footprints could not be reconstructed.

if nargin < 3 || isempty(params)
    params = struct('peri_win', [-2 3], 'z_thresh', 2);
end
lm = fullfile(out_root, 'lab_meeting');
per = fullfile(lm, 'per_recording');
avg = fullfile(lm, 'cohort_average');
if exist(per, 'dir') ~= 7; mkdir(per); end
if exist(avg, 'dir') ~= 7; mkdir(avg); end

n = numel(results);
fid = fopen(fullfile(lm, 'recording_index.csv'), 'w');
fprintf(fid, 'label,movie,mouse,session,power_mWcm2,mirror,n_neurons,n_pokes,frac_exc,frac_inh,p_pop,src_folder\n');

for i = 1:n
    R = results{i};
    info = parse_rec_label(R.movie);
    src = fullfile(out_root, R.movie);
    dst = fullfile(per, info.label);
    if exist(dst, 'dir') ~= 7; mkdir(dst); end
    copy_if(fullfile(src, 'fig1_stim_population.png'),     fullfile(dst, 'fig1_population_psth.png'));
    copy_if(fullfile(src, 'fig2_neuron_sta_heatmap.png'),  fullfile(dst, 'fig2_neuron_sta.png'));
    copy_if(fullfile(src, 'fig3_responsive_cells.png'),    fullfile(dst, 'fig3_responsive_cells.png'));
    src4 = fullfile(src, 'fig4_spatial_response.png');
    if exist(src4, 'file') ~= 2; src4 = fullfile(src, 'fig4_summary_image.png'); end
    copy_if(src4, fullfile(dst, 'fig4_spatial_response.png'));
    if contains(src4, 'summary_image')
        copy_if(src4, fullfile(dst, 'fig4_fov_anatomy_NOT_response_map.png'));
    end
    copy_if(fullfile(src, 'fig5_trial_consistency.png'),   fullfile(dst, 'fig5_trial_timing.png'));
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%d,%d,%.4f,%.4f,%.3e,%s\n', ...
        info.label, R.movie, info.mouse, info.session, info.power, info.mirror, ...
        R.n_neurons, R.n_pokes_total, R.frac_exc, R.frac_inh, R.p_pop, src);
end
fclose(fid);

write_cohort_five(results, avg, params);
write_readme(lm);
fprintf('Lab-meeting pack -> %s\n', lm);
end

function copy_if(src, dst)
if exist(src, 'file') == 2
    copyfile(src, dst);
end
end

function info = parse_rec_label(movie)
%PARSE_REC_LABEL  Friendly names + power for known 20250813 sessions.
info.mouse = 'unknown'; info.session = ''; info.mirror = ''; info.power = '';
if contains(movie, 'mirror1'); info.mirror = 'mirror1'; else; info.mirror = 'mirror2'; end
tok = regexp(movie, 's(\d+)', 'tokens', 'once');
if ~isempty(tok); info.session = ['s' tok{1}]; end

if contains(movie, '_m2c-') || contains(movie, 'm201_')
    info.mouse = 'm2c';
elseif contains(movie, 'm3_')
    info.mouse = 'm3_m1f';
elseif contains(movie, 'm200_') || contains(movie, '_m1f-') || contains(movie, '_m1_')
    info.mouse = 'm1f';
end

% 20250813 powers from Peiyi/user. Do not assume the same mW for other days.
if contains(movie, 'd250813')
    if strcmp(info.session, 's00'); info.power = '8.38';
    elseif strcmp(info.session, 's01'); info.power = '20';
    else; info.power = 'unknown';
    end
else
    info.power = 'unknown';
end

pwr = info.power;
if strcmp(pwr, 'unknown'); pwr = 'pwrNA'; else; pwr = ['pwr' pwr]; end
info.label = sprintf('%s_%s_%s_%s', info.mouse, info.session, pwr, info.mirror);
end

function write_cohort_five(results, avg, params)
n = numel(results);
tq = params.peri_win(1) : 0.05 : params.peri_win(2);
psth = nan(n, numel(tq)); frac = nan(n, numel(tq));
staA = nan(n, numel(tq)); staE = nan(n, numel(tq)); staI = nan(n, numel(tq));
col_exc = [0.85 0.20 0.15]; col_inh = [0.10 0.35 0.85];

for i = 1:n
    R = results{i};
    psth(i, :) = interp1(R.t_peri, R.pop_psth, tq, 'linear', NaN);
    frac(i, :) = interp1(R.t_peri, R.frac_psth, tq, 'linear', NaN);
    if isfield(R, 'sta_all'); staA(i, :) = interp1(R.t_peri, R.sta_all, tq, 'linear', NaN); end
    if isfield(R, 'sta_exc') && all(isfinite(R.sta_exc)); staE(i, :) = interp1(R.t_peri, R.sta_exc, tq, 'linear', NaN); end
    if isfield(R, 'sta_inh') && all(isfinite(R.sta_inh)); staI(i, :) = interp1(R.t_peri, R.sta_inh, tq, 'linear', NaN); end
end
sem = @(x) std(x, 0, 1, 'omitnan') ./ sqrt(max(sum(isfinite(x), 1), 1));

% ----- fig1 cohort population -----
figure('Position', [30 30 1500 560], 'Color', 'w');
sgtitle(sprintf('Fig1  |  All recordings (n=%d): population peri-poke', n), ...
    'FontSize', 13, 'FontWeight', 'bold');
subplot(1, 2, 1);
plot(tq, psth', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7); hold on;
mu = mean(psth, 1, 'omitnan'); se = sem(psth);
fill([tq fliplr(tq)], [mu-se fliplr(mu+se)], col_exc, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(tq, mu, 'r', 'LineWidth', 2.5); xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel('Population mean z');
title('Gray = each movie  |  red = mean \pm SEM'); box off;
subplot(1, 2, 2);
plot(tq, frac', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7); hold on;
mu = mean(frac, 1, 'omitnan'); se = sem(frac);
fill([tq fliplr(tq)], [mu-se fliplr(mu+se)], [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(tq, mu, 'k', 'LineWidth', 2.5); xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel(sprintf('Fraction active (z>%g)', params.z_thresh));
title('Population recruitment'); box off;
exportgraphics(gcf, fullfile(avg, 'fig1_population_psth.png'), 'Resolution', 200);

% ----- fig2 mean STA -----
figure('Position', [30 30 900 520], 'Color', 'w');
plot(tq, mean(staE, 1, 'omitnan'), 'Color', col_exc, 'LineWidth', 2.2); hold on;
plot(tq, mean(staI, 1, 'omitnan'), 'Color', col_inh, 'LineWidth', 2.2);
plot(tq, mean(staA, 1, 'omitnan'), 'k', 'LineWidth', 1.8);
xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel('Mean STA (z)');
title(sprintf('Fig2  |  Cohort-mean STA (average of %d recordings)', n));
legend('Excited cells (mean of movies)', 'Inhibited cells', 'All cells', 'Location', 'best');
box off;
exportgraphics(gcf, fullfile(avg, 'fig2_neuron_sta.png'), 'Resolution', 200);

% ----- fig3 responsive fractions -----
fe = nan(n, 1); fi = nan(n, 1); pp = nan(n, 1); npk = nan(n, 1); nn = nan(n, 1);
for i = 1:n
    fe(i) = results{i}.frac_exc; fi(i) = results{i}.frac_inh;
    pp(i) = results{i}.p_pop; npk(i) = results{i}.n_pokes_total; nn(i) = results{i}.n_neurons;
end
figure('Position', [30 30 1400 500], 'Color', 'w');
sgtitle('Fig3  |  Responsive cells across recordings', 'FontSize', 13, 'FontWeight', 'bold');
subplot(1, 3, 1);
b = bar([1 2], [mean(fe) mean(fi)]*100, 0.55, 'EdgeColor', 'k'); hold on;
b.FaceColor = 'flat'; b.CData = [col_exc; col_inh];
errorbar([1 2], [mean(fe) mean(fi)]*100, [std(fe) std(fi)]*100, 'k.', 'LineWidth', 1.4);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Excited', 'Inhibited'});
ylabel('% of neurons (mean \pm SD across movies)'); box off;
title(sprintf('n = %d recordings', n));
subplot(1, 3, 2);
plot(100*fe, 100*fi, 'ko', 'MarkerFaceColor', [0.3 0.3 0.3]); hold on;
xlabel('% excited'); ylabel('% inhibited'); axis equal;
title('Each movie'); box off; grid on;
subplot(1, 3, 3);
histogram(log10(max(pp, realmin)), 'FaceColor', [0.5 0.5 0.7], 'EdgeColor', 'none');
xlabel('log10 population p'); ylabel('# movies');
title(sprintf('%d / %d movies p<0.05', sum(pp<0.05), n)); box off;
exportgraphics(gcf, fullfile(avg, 'fig3_responsive_cells.png'), 'Resolution', 200);

% ----- fig4 spatial response gallery -----
has_maps = false;
for i = 1:n
    if isfield(results{i}, 'resp_map_small') && ~isempty(results{i}.resp_map_small)
        has_maps = true;
        break
    end
end
figure('Position', [20 20 1600 900], 'Color', 'w');
if has_maps
    sgtitle({'Fig4  |  Poke response on FOV  (red = excited, blue = inhibited)', ...
        'Each panel is one movie; FOVs are not spatially aligned across mice/mirrors'}, ...
        'FontSize', 12, 'FontWeight', 'bold');
else
    sgtitle({'Fig4  |  FOV anatomy (spatial footprints unavailable)', ...
        'Summary image only — not a poke on/off map'}, ...
        'FontSize', 12, 'FontWeight', 'bold');
end
ncol = 4; nrow = ceil(n / ncol);
half_n = 128;
cmap_rb = [[linspace(0, 1, half_n) ones(1, half_n)]' ...
           [linspace(0, 1, half_n) linspace(1, 0, half_n)]' ...
           [ones(1, half_n) linspace(1, 0, half_n)]'];
for i = 1:n
    subplot(nrow, ncol, i);
    R = results{i};
    info = parse_rec_label(R.movie);
    if has_maps && isfield(R, 'resp_map_small') && ~isempty(R.resp_map_small)
        M = R.resp_map_small;
        climv = prctile(abs(M(M ~= 0)), 98);
        if isempty(climv) || ~isfinite(climv) || climv == 0; climv = 1; end
        imagesc(M); colormap(gca, cmap_rb); clim([-climv climv]); axis image off;
    elseif isfield(R, 'sumimg_small') && ~isempty(R.sumimg_small)
        imagesc(R.sumimg_small); colormap(gca, 'gray'); axis image off;
    else
        axis off; text(0.5, 0.5, 'no image', 'HorizontalAlignment', 'center');
    end
    title(strrep(info.label, '_', '\_'), 'FontSize', 8, 'Interpreter', 'tex');
end
exportgraphics(gcf, fullfile(avg, 'fig4_spatial_response.png'), 'Resolution', 180);

% ----- fig5 timing -----
alli = [];
for i = 1:n
    if isfield(results{i}, 'iti') && ~isempty(results{i}.iti)
        alli = [alli; results{i}.iti(:)]; %#ok<AGROW>
    end
end
figure('Position', [30 30 1400 500], 'Color', 'w');
sgtitle('Fig5  |  Stimulus counts and inter-poke interval', 'FontSize', 13, 'FontWeight', 'bold');
subplot(1, 2, 1);
bar(1:n, npk, 0.7, 'FaceColor', [0.4 0.45 0.75], 'EdgeColor', 'none');
ylabel('# detected pokes'); xlabel('Recording');
title(sprintf('Pokes per movie (mean %.1f)', mean(npk))); box off;
subplot(1, 2, 2);
if ~isempty(alli)
    histogram(alli, 'FaceColor', [0.4 0.4 0.7], 'EdgeColor', 'none');
    xlabel('Inter-poke interval (s)'); ylabel('Count');
    title(sprintf('Pooled ITI  median = %.2f s', median(alli)));
end
box off;
exportgraphics(gcf, fullfile(avg, 'fig5_trial_timing.png'), 'Resolution', 200);
end

function write_readme(lm)
fid = fopen(fullfile(lm, 'README.txt'), 'w');
fprintf(fid, ['Lab-meeting figure pack\n' ...
    '========================\n\n' ...
    'per_recording/<mouse_session_power_mirror>/\n' ...
    '  fig1 population PSTH for that movie\n' ...
    '  fig2 per-neuron STA heatmap\n' ...
    '  fig3 excited / inhibited counts\n' ...
    '  fig4 poke response on FOV: red = excited, blue = inhibited\n' ...
    '  fig5 per-trial raster + ITI\n\n' ...
    'cohort_average/\n' ...
    '  Same 5 figure numbers, averaged (or pooled) across all recordings.\n' ...
    '  fig4 is a gallery of per-movie response maps (FOVs not aligned).\n']);
fclose(fid);
end

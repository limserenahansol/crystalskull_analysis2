% EXTRACT Results: Baseline vs Injection — Population-Level Analysis
% Crystal-skull / Mesoscope whole dorsal cortex calcium imaging
%
% Run once — all figures and summary saved to out_dir (see CONFIG below).
%
% See README.md, docs/DATASET_STRUCTURE.md, and docs/METHODS.md (metrics & statistics).

clear; close all; clc;
tic;

%% ========== CONFIG — Edit these paths for your data ==========
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\CS1014-1b';
bl_subfolder  = 'mirror1path';   % subfolder under Baseline/
inj_subfolder = 'mirror1path';   % subfolder under Injection/
out_dir       = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures';

% Analysis parameters
z_thresh      = 2;
n_boot        = 1000;
time_bin      = 1;              % seconds for population binning
n_pc          = 10;
n_corr_sample = 500;
frame_px      = 2304;
max_show      = 2000;
%% =============================================================

bl_path  = fullfile(base_dir, 'Baseline',  bl_subfolder);
inj_path = fullfile(base_dir, 'Injection', inj_subfolder);

col_bl  = [0.13 0.59 0.95];
col_inj = [1.00 0.34 0.13];

% --- Load data (using h5read for HDF5/v7.3 mat files) --------------------
frr_bl_file  = fullfile(bl_path,  'M_moco_frr.mat');
frr_inj_file = fullfile(inj_path, 'M_moco_frr.mat');
ext_bl_file  = fullfile(bl_path,  'M_moco_ds_ext.mat');
ext_inj_file = fullfile(inj_path, 'M_moco_ds_ext.mat');
sum_bl_file  = fullfile(bl_path,  'M_summary.mat');
sum_inj_file = fullfile(inj_path, 'M_summary.mat');

fprintf('Loading Baseline frr (T + S)...\n');
T_bl = double(h5read(frr_bl_file, '/T'))';
S_bl_ir   = double(h5read(frr_bl_file, '/S/ir'));
S_bl_jc   = double(h5read(frr_bl_file, '/S/jc'));
S_bl_data = double(h5read(frr_bl_file, '/S/data'));

fprintf('Loading Injection frr (T + S)...\n');
T_inj = double(h5read(frr_inj_file, '/T'))';
S_inj_ir   = double(h5read(frr_inj_file, '/S/ir'));
S_inj_jc   = double(h5read(frr_inj_file, '/S/jc'));
S_inj_data = double(h5read(frr_inj_file, '/S/data'));

fprintf('Loading EXTRACT temporal weights...\n');
tw_bl  = double(h5read(ext_bl_file,  '/temporal_weights'))';
tw_inj = double(h5read(ext_inj_file, '/temporal_weights'))';
sumimg_bl  = h5read(ext_bl_file,  '/info/summary_image')';
sumimg_inj = h5read(ext_inj_file, '/info/summary_image')';

fprintf('Loading framerates...\n');
fs_bl  = h5read(sum_bl_file,  '/framerate');
fs_inj = h5read(sum_inj_file, '/framerate');

fprintf('Loading cell map images...\n');
cellmap_bl  = imread(fullfile(bl_path,  'cell_map.png'));
cellmap_inj = imread(fullfile(inj_path, 'cell_map.png'));

[n_frames_bl,  n_neurons_bl]  = size(T_bl);
[n_frames_inj, n_neurons_inj] = size(T_inj);
dur_bl  = n_frames_bl  / fs_bl;
dur_inj = n_frames_inj / fs_inj;

fprintf('\n  Baseline:  %d neurons | %d frames | %.2f Hz | %.1f s\n', n_neurons_bl, n_frames_bl, fs_bl, dur_bl);
fprintf('  Injection: %d neurons | %d frames | %.2f Hz | %.1f s\n\n', n_neurons_inj, n_frames_inj, fs_inj, dur_inj);

% --- Compute dF/F (inline, 8th percentile rolling baseline) ---------------
fprintf('Computing dF/F (Baseline)...\n');
nf = n_frames_bl; win = max(round(nf*0.15),300); half = floor(win/2);
F0 = nan(size(T_bl));
for c = 1:half:nf
    s = max(1,c-half); e = min(nf,c+half);
    f0_chunk = prctile(T_bl(s:e,:),8,1);
    fe = min(c+half-1,nf);
    F0(c:fe,:) = repmat(f0_chunk, fe-c+1, 1);
end
for j = 1:size(F0,2)
    col = F0(:,j); m = isnan(col);
    if any(m), iv = find(~m); col(m) = interp1(iv,col(iv),find(m),'nearest','extrap'); F0(:,j) = col; end
end
F0(F0<0.5)=0.5;
dff_bl = (T_bl - F0) ./ F0;
clear F0

fprintf('Computing dF/F (Injection)...\n');
nf = n_frames_inj; win = max(round(nf*0.15),300); half = floor(win/2);
F0 = nan(size(T_inj));
for c = 1:half:nf
    s = max(1,c-half); e = min(nf,c+half);
    f0_chunk = prctile(T_inj(s:e,:),8,1);
    fe = min(c+half-1,nf);
    F0(c:fe,:) = repmat(f0_chunk, fe-c+1, 1);
end
for j = 1:size(F0,2)
    col = F0(:,j); m = isnan(col);
    if any(m), iv = find(~m); col(m) = interp1(iv,col(iv),find(m),'nearest','extrap'); F0(:,j) = col; end
end
F0(F0<0.5)=0.5;
dff_inj = (T_inj - F0) ./ F0;
clear F0 T_bl T_inj

% --- Z-score each neuron -------------------------------------------------
fprintf('Z-scoring...\n');
z_bl  = (dff_bl  - mean(dff_bl,1))  ./ max(std(dff_bl,0,1),  1e-8);
z_inj = (dff_inj - mean(dff_inj,1)) ./ max(std(dff_inj,0,1), 1e-8);

% --- Per-neuron statistics (inline) ---------------------------------------
fprintf('Per-neuron statistics...\n');
mean_dff_bl   = mean(dff_bl,1)';      mean_dff_inj   = mean(dff_inj,1)';
std_dff_bl    = std(dff_bl,0,1)';     std_dff_inj    = std(dff_inj,0,1)';
max_dff_bl    = max(dff_bl,[],1)';    max_dff_inj    = max(dff_inj,[],1)';
skew_dff_bl   = skewness(dff_bl,1,1)'; skew_dff_inj  = skewness(dff_inj,1,1)';

cross_bl  = diff(double(z_bl  > z_thresh),1,1);
cross_inj = diff(double(z_inj > z_thresh),1,1);
event_count_bl  = sum(cross_bl  == 1, 1)';
event_count_inj = sum(cross_inj == 1, 1)';
event_rate_bl   = event_count_bl  / dur_bl;
event_rate_inj  = event_count_inj / dur_inj;
clear cross_bl cross_inj

% --- Population vectors (1-s bins) ----------------------------------------
fprintf('Population rate vectors...\n');
bin_frames_bl  = round(fs_bl  * time_bin);
bin_frames_inj = round(fs_inj * time_bin);
n_bins_bl  = floor(n_frames_bl  / bin_frames_bl);
n_bins_inj = floor(n_frames_inj / bin_frames_inj);

pop_rate_bl  = zeros(n_bins_bl,1);  pop_frac_bl  = zeros(n_bins_bl,1);
pop_rate_inj = zeros(n_bins_inj,1); pop_frac_inj = zeros(n_bins_inj,1);
t_pop_bl  = ((1:n_bins_bl)'-0.5)  * time_bin;
t_pop_inj = ((1:n_bins_inj)'-0.5) * time_bin;

for b = 1:n_bins_bl
    idx = (b-1)*bin_frames_bl+1 : b*bin_frames_bl;
    chunk = z_bl(idx,:);
    pop_rate_bl(b) = mean(chunk(:));
    pop_frac_bl(b) = mean(any(chunk > z_thresh, 1));
end
for b = 1:n_bins_inj
    idx = (b-1)*bin_frames_inj+1 : b*bin_frames_inj;
    chunk = z_inj(idx,:);
    pop_rate_inj(b) = mean(chunk(:));
    pop_frac_inj(b) = mean(any(chunk > z_thresh, 1));
end

% Wilcoxon: compare per-1s-bin fraction-active distributions (fig2 synchrony panel + summary_stats)
[p_popfrac,~] = ranksum(pop_frac_bl, pop_frac_inj);

% --- Pairwise correlations (sampled) --------------------------------------
fprintf('Pairwise correlations (%d sampled)...\n', n_corr_sample);
rng(42);
idx_c_bl  = randperm(n_neurons_bl,  min(n_corr_sample, n_neurons_bl));
idx_c_inj = randperm(n_neurons_inj, min(n_corr_sample, n_neurons_inj));
C_bl  = corrcoef(z_bl(:, idx_c_bl));
C_inj = corrcoef(z_inj(:, idx_c_inj));
mask_triu = triu(true(size(C_bl,1)), 1);
pw_bl  = C_bl(mask_triu);
pw_inj = C_inj(mask_triu);

% --- PCA on binned population activity ------------------------------------
fprintf('PCA on population activity...\n');
Z_bin_bl  = zeros(n_bins_bl,  n_neurons_bl);
Z_bin_inj = zeros(n_bins_inj, n_neurons_inj);
for b = 1:n_bins_bl
    Z_bin_bl(b,:) = mean(z_bl((b-1)*bin_frames_bl+1 : b*bin_frames_bl, :), 1);
end
for b = 1:n_bins_inj
    Z_bin_inj(b,:) = mean(z_inj((b-1)*bin_frames_inj+1 : b*bin_frames_inj, :), 1);
end
[~, score_bl,  ~, ~, expl_bl]  = pca(Z_bin_bl,  'NumComponents', n_pc);
[~, score_inj, ~, ~, expl_inj] = pca(Z_bin_inj, 'NumComponents', n_pc);
clear Z_bin_bl Z_bin_inj

% --- Spatial activity maps (project metric onto FOV via sparse S) ----------
fprintf('Spatial activity maps...\n');
n_pix = frame_px * frame_px;

% Baseline maps
jc = S_bl_jc; col_idx = zeros(jc(end),1);
for c = 1:n_neurons_bl, col_idx(jc(c)+1:jc(c+1)) = c; end
Ss = sparse(S_bl_ir+1, col_idx, S_bl_data, n_pix, n_neurons_bl);
cs = full(sum(Ss,1)); cs(cs==0)=1; Sn = Ss ./ cs;
actmap_bl     = reshape(full(Sn * mean_dff_bl),   frame_px, frame_px);
erate_map_bl  = reshape(full(Sn * event_rate_bl),  frame_px, frame_px);
clear Ss Sn col_idx jc cs S_bl_ir S_bl_jc S_bl_data

% Injection maps
jc = S_inj_jc; col_idx = zeros(jc(end),1);
for c = 1:n_neurons_inj, col_idx(jc(c)+1:jc(c+1)) = c; end
Ss = sparse(S_inj_ir+1, col_idx, S_inj_data, n_pix, n_neurons_inj);
cs = full(sum(Ss,1)); cs(cs==0)=1; Sn = Ss ./ cs;
actmap_inj    = reshape(full(Sn * mean_dff_inj),   frame_px, frame_px);
erate_map_inj = reshape(full(Sn * event_rate_inj), frame_px, frame_px);
clear Ss Sn col_idx jc cs S_inj_ir S_inj_jc S_inj_data

clear dff_bl dff_inj

% --- Build colormaps ------------------------------------------------------
t_cm = linspace(0,1,256)';
inf_r = min(1,max(0, -0.0155+t_cm.*(5.27+t_cm.*(-14.35+t_cm.*(17.77+t_cm.*(-9.42+t_cm.*1.76))))));
inf_g = min(1,max(0, -0.0123+t_cm.*(-0.23+t_cm.*(7.77+t_cm.*(-17.05+t_cm.*(16.95+t_cm.*(-6.34)))))));
inf_b = min(1,max(0,  0.095 +t_cm.*(4.00+t_cm.*(-15.24+t_cm.*(23.56+t_cm.*(-17.39+t_cm.*5.01))))));
cmap_inf = [inf_r inf_g inf_b];

half_n = 128;
cmap_rb = [[linspace(0,1,half_n) ones(1,half_n)]' ...
           [linspace(0,1,half_n) linspace(1,0,half_n)]' ...
           [ones(1,half_n) linspace(1,0,half_n)]'];

t_v = linspace(0,1,256)';
vr = min(1,max(0, 0.267+t_v.*(-0.005+t_v.*(3.07+t_v.*(-6.36+t_v.*(5.27+t_v.*(-1.58)))))));
vg = min(1,max(0, 0.004+t_v.*(1.42+t_v.*(-0.15+t_v.*(-1.95+t_v.*(3.11+t_v.*(-1.40)))))));
vb = min(1,max(0, 0.329+t_v.*(1.44+t_v.*(-5.10+t_v.*(9.18+t_v.*(-8.28+t_v.*2.75))))));
cmap_viridis = [vr vg vb];

fprintf('\nPlotting figures...\n');

% --- FIGURE 1: Cell Maps & Summary Images ---------------------------------
figure('Position',[30 30 1700 900],'Color','w');
sgtitle('Field of View & Detected Neurons','FontSize',15,'FontWeight','bold');

subplot(2,3,1);
imagesc(sumimg_bl); colormap(gca,'gray'); axis image off;
title(sprintf('Baseline Summary Image\n(%dx%d)', frame_px, frame_px));

subplot(2,3,4);
imagesc(sumimg_inj); colormap(gca,'gray'); axis image off;
title(sprintf('Injection Summary Image\n(%dx%d)', frame_px, frame_px));

subplot(2,3,2);
imshow(cellmap_bl);
title(sprintf('Baseline Cell Map\n%d neurons', n_neurons_bl));

subplot(2,3,5);
imshow(cellmap_inj);
title(sprintf('Injection Cell Map\n%d neurons', n_neurons_inj));

subplot(2,3,3);
imagesc(actmap_bl); colormap(gca, cmap_viridis); axis image off;
cb = colorbar; cb.Label.String = 'Mean dF/F';
title('Baseline: Spatial Mean Activity');

subplot(2,3,6);
imagesc(actmap_inj); colormap(gca, cmap_viridis); axis image off;
cb = colorbar; cb.Label.String = 'Mean dF/F';
title('Injection: Spatial Mean Activity');

exportgraphics(gcf, fullfile(out_dir,'fig1_cellmaps_FOV.png'),'Resolution',200);

% --- FIGURE 2: Population Activity Overview --------------------------------
figure('Position',[30 30 1600 1000],'Color','w');
sgtitle('Population-Level Activity: Baseline vs Injection','FontSize',15,'FontWeight','bold');

subplot(3,2,1);
plot(t_pop_bl,  pop_rate_bl,  'Color',[col_bl 0.8],  'LineWidth',0.6); hold on;
plot(t_pop_inj, pop_rate_inj, 'Color',[col_inj 0.8], 'LineWidth',0.6);
ylabel('Pop. Mean z-score'); xlabel('Time (s)');
title('Population Mean Activity Over Time');
legend('Baseline','Injection','Location','northeast'); box off;

subplot(3,2,2);
plot(t_pop_bl,  pop_frac_bl,  'Color',[col_bl 0.8],  'LineWidth',0.6); hold on;
plot(t_pop_inj, pop_frac_inj, 'Color',[col_inj 0.8], 'LineWidth',0.6);
ylabel(sprintf('Fraction Active (z>%d)', z_thresh)); xlabel('Time (s)');
title('Population Recruitment Over Time');
legend('Baseline','Injection','Location','northeast'); box off;

subplot(3,2,3);
edges_sync = linspace(0, max(max(pop_frac_bl),max(pop_frac_inj)), 50);
histogram(pop_frac_bl,  edges_sync,'Normalization','pdf','FaceColor',col_bl, 'FaceAlpha',0.5,'EdgeColor','none'); hold on;
histogram(pop_frac_inj, edges_sync,'Normalization','pdf','FaceColor',col_inj,'FaceAlpha',0.5,'EdgeColor','none');
rng(1);
boot_bl  = bootstrp(n_boot, @mean, pop_frac_bl);
boot_inj = bootstrp(n_boot, @mean, pop_frac_inj);
ci_bl  = prctile(boot_bl,  [2.5 97.5]);
ci_inj = prctile(boot_inj, [2.5 97.5]);
xline(mean(pop_frac_bl),  '-','Color',col_bl, 'LineWidth',2);
xline(mean(pop_frac_inj), '-','Color',col_inj,'LineWidth',2);
xlabel('Fraction Active'); ylabel('Density');
title(sprintf(['Synchrony Distribution\nBL: %.4f [%.4f-%.4f]  Inj: %.4f [%.4f-%.4f]\n' ...
    'Wilcoxon rank-sum (per-bin fractions), p = %.2e'], ...
    mean(pop_frac_bl), ci_bl(1), ci_bl(2), mean(pop_frac_inj), ci_inj(1), ci_inj(2), p_popfrac));
legend('Baseline','Injection'); box off;

subplot(3,2,4);
rng(2);
boot_er_bl  = bootstrp(n_boot, @mean, event_rate_bl);
boot_er_inj = bootstrp(n_boot, @mean, event_rate_inj);
ci_er_bl  = prctile(boot_er_bl,  [2.5 97.5]);
ci_er_inj = prctile(boot_er_inj, [2.5 97.5]);
bar_x = [1 2];
bar_y = [mean(event_rate_bl) mean(event_rate_inj)];
bar_err = [bar_y(1)-ci_er_bl(1), bar_y(2)-ci_er_inj(1); ci_er_bl(2)-bar_y(1), ci_er_inj(2)-bar_y(2)];
b = bar(bar_x, bar_y, 0.5, 'EdgeColor','k'); hold on;
b.FaceColor = 'flat'; b.CData = [col_bl; col_inj];
errorbar(bar_x, bar_y, bar_err(1,:), bar_err(2,:), 'k.', 'LineWidth',1.5, 'CapSize',10);
set(gca,'XTick',[1 2],'XTickLabel',{'Baseline','Injection'});
ylabel('Mean Event Rate (events/s)');
[p_er,~] = ranksum(event_rate_bl, event_rate_inj);
title(sprintf('Population Mean Event Rate\n95%% CI, p = %.2e', p_er));
yl = ylim; ylim([0 yl(2)*1.15]);
line([1 2],[yl(2)*1.02 yl(2)*1.02],'Color','k','LineWidth',1.5);
text(1.5, yl(2)*1.07, sprintf('***  p = %.1e', p_er), 'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
box off;

subplot(3,2,5);
rng(3);
boot_dff_bl  = bootstrp(n_boot, @mean, mean_dff_bl);
boot_dff_inj = bootstrp(n_boot, @mean, mean_dff_inj);
ci_dff_bl  = prctile(boot_dff_bl,  [2.5 97.5]);
ci_dff_inj = prctile(boot_dff_inj, [2.5 97.5]);
bar_y2 = [mean(mean_dff_bl) mean(mean_dff_inj)];
bar_err2 = [bar_y2(1)-ci_dff_bl(1), bar_y2(2)-ci_dff_inj(1); ci_dff_bl(2)-bar_y2(1), ci_dff_inj(2)-bar_y2(2)];
b2 = bar(bar_x, bar_y2, 0.5,'EdgeColor','k'); hold on;
b2.FaceColor = 'flat'; b2.CData = [col_bl; col_inj];
errorbar(bar_x, bar_y2, bar_err2(1,:), bar_err2(2,:), 'k.','LineWidth',1.5,'CapSize',10);
set(gca,'XTick',[1 2],'XTickLabel',{'Baseline','Injection'});
ylabel('Mean dF/F');
[p_dff,~] = ranksum(mean_dff_bl, mean_dff_inj);
title(sprintf('Population Mean dF/F\n95%% CI, p = %.2e', p_dff));
yl = ylim; ylim([0 yl(2)*1.15]);
line([1 2],[yl(2)*1.02 yl(2)*1.02],'Color','k','LineWidth',1.5);
if p_dff < 0.001, sig_str = '***'; elseif p_dff < 0.01, sig_str = '**'; elseif p_dff < 0.05, sig_str = '*'; else, sig_str = 'n.s.'; end
text(1.5, yl(2)*1.07, sprintf('%s  p = %.1e', sig_str, p_dff), 'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
box off;

subplot(3,2,6);
b3 = bar([1 2], [n_neurons_bl n_neurons_inj], 0.5,'EdgeColor','k');
b3.FaceColor = 'flat'; b3.CData = [col_bl; col_inj];
set(gca,'XTick',[1 2],'XTickLabel',{'Baseline','Injection'});
ylabel('Neurons'); title('Detected Neurons');
text(1, n_neurons_bl+200,  sprintf('%d', n_neurons_bl),  'HorizontalAlignment','center','FontWeight','bold');
text(2, n_neurons_inj+200, sprintf('%d', n_neurons_inj), 'HorizontalAlignment','center','FontWeight','bold');
box off;

exportgraphics(gcf, fullfile(out_dir,'fig2_population_activity.png'),'Resolution',200);

% --- FIGURE 3: Per-Neuron CDF Distributions --------------------------------
figure('Position',[30 30 1600 500],'Color','w');
sgtitle('Per-Neuron Metric Distributions (CDF)','FontSize',14,'FontWeight','bold');

cdf_bl  = {mean_dff_bl,  event_rate_bl,  std_dff_bl,  skew_dff_bl};
cdf_inj = {mean_dff_inj, event_rate_inj, std_dff_inj, skew_dff_inj};
cdf_labels = {'Mean dF/F', 'Event Rate (events/s)', 'Std dF/F', 'Skewness'};

for i = 1:4
    subplot(1,4,i);
    d1 = cdf_bl{i}; d2 = cdf_inj{i};
    [f1,x1] = ecdf(d1); [f2,x2] = ecdf(d2);
    plot(x1,f1,'Color',col_bl,'LineWidth',2); hold on;
    plot(x2,f2,'Color',col_inj,'LineWidth',2);
    xlabel(cdf_labels{i}); ylabel('CDF');
    [~, p_ks] = kstest2(d1, d2);
    [p_rs,~]  = ranksum(d1, d2);
    title(sprintf('%s\nKS p=%.1e | RS p=%.1e', cdf_labels{i}, p_ks, p_rs), 'FontSize',9);
    xl = prctile([d1; d2], [1 99]); xlim(xl);
    legend(sprintf('BL (med=%.3f)', median(d1)), ...
           sprintf('Inj (med=%.3f)', median(d2)), 'Location','southeast','FontSize',7);
    box off;
end

exportgraphics(gcf, fullfile(out_dir,'fig3_neuron_distributions.png'),'Resolution',200);

% --- FIGURE 4: Pairwise Correlations ---------------------------------------
figure('Position',[30 30 1600 500],'Color','w');
sgtitle(sprintf('Pairwise Neuron Correlations (n=%d sampled)', n_corr_sample),'FontSize',14,'FontWeight','bold');

subplot(1,3,1);
imagesc(C_bl,[-0.3 0.3]); colormap(gca, cmap_rb);
colorbar; axis image; xlabel('Neuron'); ylabel('Neuron');
title('Baseline');

subplot(1,3,2);
imagesc(C_inj,[-0.3 0.3]); colormap(gca, cmap_rb);
colorbar; axis image; xlabel('Neuron');
title('Injection');

subplot(1,3,3);
edges_c = linspace(-0.3, 0.5, 80);
histogram(pw_bl,  edges_c,'Normalization','pdf','FaceColor',col_bl, 'FaceAlpha',0.5,'EdgeColor','none'); hold on;
histogram(pw_inj, edges_c,'Normalization','pdf','FaceColor',col_inj,'FaceAlpha',0.5,'EdgeColor','none');
xline(median(pw_bl), '--','Color',col_bl, 'LineWidth',2);
xline(median(pw_inj),'--','Color',col_inj,'LineWidth',2);
xlabel('Pearson r'); ylabel('Density');
[p_pw,~] = ranksum(pw_bl, pw_inj);
title(sprintf('Pairwise Correlation\nBL med=%.4f  Inj med=%.4f  p=%.2e', median(pw_bl), median(pw_inj), p_pw));
legend('Baseline','Injection','Location','northeast'); box off;

exportgraphics(gcf, fullfile(out_dir,'fig4_network_correlations.png'),'Resolution',200);

% --- FIGURE 5: PCA / Dimensionality ----------------------------------------
figure('Position',[30 30 1600 500],'Color','w');
sgtitle('Population Dimensionality (PCA)','FontSize',14,'FontWeight','bold');

subplot(1,3,1);
plot(1:n_pc, cumsum(expl_bl(1:n_pc)),  'o-','Color',col_bl, 'LineWidth',2,'MarkerFaceColor',col_bl); hold on;
plot(1:n_pc, cumsum(expl_inj(1:n_pc)), 's-','Color',col_inj,'LineWidth',2,'MarkerFaceColor',col_inj);
xlabel('# PCs'); ylabel('Cumulative % Variance');
title('Scree Plot (cumulative)');
legend(sprintf('BL (PC1=%.1f%%)', expl_bl(1)), sprintf('Inj (PC1=%.1f%%)', expl_inj(1)), 'Location','southeast');
box off; grid on;

subplot(1,3,2);
scatter(score_bl(:,1), score_bl(:,2), 8, linspace(0,1,n_bins_bl)', 'filled');
colormap(gca, cmap_inf); cb = colorbar; cb.Label.String = 'Normalized Time';
xlabel(sprintf('PC1 (%.1f%%)', expl_bl(1)));
ylabel(sprintf('PC2 (%.1f%%)', expl_bl(2)));
title('Baseline: Neural Trajectory'); box off;

subplot(1,3,3);
scatter(score_inj(:,1), score_inj(:,2), 8, linspace(0,1,n_bins_inj)', 'filled');
colormap(gca, cmap_inf); cb = colorbar; cb.Label.String = 'Normalized Time';
xlabel(sprintf('PC1 (%.1f%%)', expl_inj(1)));
ylabel(sprintf('PC2 (%.1f%%)', expl_inj(2)));
title('Injection: Neural Trajectory'); box off;

exportgraphics(gcf, fullfile(out_dir,'fig5_PCA_dimensionality.png'),'Resolution',200);

% --- FIGURE 6: Activity Heatmaps ------------------------------------------
figure('Position',[30 30 1600 650],'Color','w');
sgtitle('Z-scored Activity Heatmaps (2000 sampled neurons, sorted by peak)','FontSize',14,'FontWeight','bold');

subplot(1,2,1);
rng(99);
if n_neurons_bl > max_show, idx_h = randperm(n_neurons_bl, max_show); else, idx_h = 1:n_neurons_bl; end
z_sub = z_bl(:, idx_h);
[~,pk] = max(z_sub,[],1); [~,ord] = sort(pk);
imagesc([0 dur_bl],[1 numel(idx_h)], z_sub(:,ord)');
caxis([0 5]); colormap(gca, cmap_inf);
xlabel('Time (s)'); ylabel('Neuron');
title(sprintf('Baseline (%d neurons)', n_neurons_bl));

subplot(1,2,2);
rng(99);
if n_neurons_inj > max_show, idx_h = randperm(n_neurons_inj, max_show); else, idx_h = 1:n_neurons_inj; end
z_sub = z_inj(:, idx_h);
[~,pk] = max(z_sub,[],1); [~,ord] = sort(pk);
imagesc([0 dur_inj],[1 numel(idx_h)], z_sub(:,ord)');
caxis([0 5]); colormap(gca, cmap_inf);
cb = colorbar; cb.Label.String = 'Z-score';
xlabel('Time (s)'); ylabel('Neuron');
title(sprintf('Injection (%d neurons)', n_neurons_inj));

exportgraphics(gcf, fullfile(out_dir,'fig6_activity_heatmaps.png'),'Resolution',200);

% --- FIGURE 7: Spatial Event Rate Maps -------------------------------------
figure('Position',[30 30 1400 600],'Color','w');
sgtitle('Spatial Event Rate Maps (events/s projected onto FOV)','FontSize',14,'FontWeight','bold');

clim_er = [0 prctile([erate_map_bl(erate_map_bl>0); erate_map_inj(erate_map_inj>0)], 95)];

subplot(1,2,1);
imagesc(erate_map_bl); colormap(gca, cmap_viridis);
caxis(clim_er); axis image off;
cb = colorbar; cb.Label.String = 'events/s';
title('Baseline: Event Rate');

subplot(1,2,2);
imagesc(erate_map_inj); colormap(gca, cmap_viridis);
caxis(clim_er); axis image off;
cb = colorbar; cb.Label.String = 'events/s';
title('Injection: Event Rate');

exportgraphics(gcf, fullfile(out_dir,'fig7_spatial_event_rate.png'),'Resolution',200);

% --- Summary statistics ----------------------------------------------------
fprintf('\nSUMMARY STATISTICS\n');
fprintf('%s\n', repmat('=', 1, 90));

fid = fopen(fullfile(out_dir, 'summary_stats.txt'), 'w');
L = {};
L{end+1} = 'EXTRACT Results: Baseline vs Injection';
L{end+1} = repmat('=', 1, 90);
L{end+1} = sprintf('%-45s %15s %15s', 'Metric', 'Baseline', 'Injection');
L{end+1} = repmat('-', 1, 75);
L{end+1} = sprintf('%-45s %15d %15d', 'Neurons detected', n_neurons_bl, n_neurons_inj);
L{end+1} = sprintf('%-45s %15d %15d', 'Frames', n_frames_bl, n_frames_inj);
L{end+1} = sprintf('%-45s %15.2f %15.2f', 'Frame rate (Hz)', fs_bl, fs_inj);
L{end+1} = sprintf('%-45s %15.1f %15.1f', 'Duration (s)', dur_bl, dur_inj);
L{end+1} = '';
L{end+1} = '--- Population-level (bootstrap 95% CI) ---';
L{end+1} = sprintf('%-45s %7.4f [%7.4f-%7.4f]    %7.4f [%7.4f-%7.4f]', ...
    'Mean event rate (events/s)', mean(event_rate_bl), ci_er_bl, mean(event_rate_inj), ci_er_inj);
L{end+1} = sprintf('%-45s %7.4f [%7.4f-%7.4f]    %7.4f [%7.4f-%7.4f]', ...
    'Mean fraction active (z>2)', mean(pop_frac_bl), ci_bl, mean(pop_frac_inj), ci_inj);
L{end+1} = sprintf('%-45s %7.4f [%7.4f-%7.4f]    %7.4f [%7.4f-%7.4f]', ...
    'Mean dF/F', mean(mean_dff_bl), ci_dff_bl, mean(mean_dff_inj), ci_dff_inj);
L{end+1} = '';
L{end+1} = '--- Per-neuron medians ---';
L{end+1} = sprintf('%-45s %15.4f %15.4f', 'Median mean dF/F',    median(mean_dff_bl),   median(mean_dff_inj));
L{end+1} = sprintf('%-45s %15.4f %15.4f', 'Median std dF/F',     median(std_dff_bl),    median(std_dff_inj));
L{end+1} = sprintf('%-45s %15.4f %15.4f', 'Median event rate',   median(event_rate_bl), median(event_rate_inj));
L{end+1} = sprintf('%-45s %15.4f %15.4f', 'Median skewness',     median(skew_dff_bl),   median(skew_dff_inj));
L{end+1} = sprintf('%-45s %15.4f %15.4f', 'Median pairwise corr', median(pw_bl), median(pw_inj));
L{end+1} = '';
L{end+1} = '--- Dimensionality (PCA) ---';
L{end+1} = sprintf('%-45s %15.2f %15.2f', 'Var explained by PC1 (%)', expl_bl(1), expl_inj(1));
L{end+1} = sprintf('%-45s %15.2f %15.2f', 'Var explained by PC1-3 (%)', sum(expl_bl(1:3)), sum(expl_inj(1:3)));
L{end+1} = '';
L{end+1} = '--- Statistical Tests ---';

[p1,~] = ranksum(mean_dff_bl,   mean_dff_inj);
[p2,~] = ranksum(event_rate_bl, event_rate_inj);
[p3,~] = ranksum(pw_bl, pw_inj);
p4 = p_popfrac;
[~,p5] = kstest2(event_rate_bl, event_rate_inj);
[~,p6] = kstest2(mean_dff_bl,   mean_dff_inj);

% Cohen's d
nx = n_neurons_bl; ny = n_neurons_inj;
sp_er  = sqrt(((nx-1)*var(event_rate_bl) + (ny-1)*var(event_rate_inj)) / (nx+ny-2));
d_er   = (mean(event_rate_bl) - mean(event_rate_inj)) / max(sp_er, 1e-12);
sp_dff = sqrt(((nx-1)*var(mean_dff_bl) + (ny-1)*var(mean_dff_inj)) / (nx+ny-2));
d_dff  = (mean(mean_dff_bl) - mean(mean_dff_inj)) / max(sp_dff, 1e-12);
np_bl = numel(pw_bl); np_inj = numel(pw_inj);
sp_pw = sqrt(((np_bl-1)*var(pw_bl) + (np_inj-1)*var(pw_inj)) / (np_bl+np_inj-2));
d_pw  = (mean(pw_bl) - mean(pw_inj)) / max(sp_pw, 1e-12);

L{end+1} = sprintf('  %-40s  p = %.2e  (Cohen d = %.3f)', 'Mean dF/F (Wilcoxon)',      p1, d_dff);
L{end+1} = sprintf('  %-40s  p = %.2e  (Cohen d = %.3f)', 'Event rate (Wilcoxon)',     p2, d_er);
L{end+1} = sprintf('  %-40s  p = %.2e  (Cohen d = %.3f)', 'Pairwise corr (Wilcoxon)',  p3, d_pw);
L{end+1} = sprintf('  %-40s  p = %.2e', 'Fraction active (Wilcoxon)', p4);
L{end+1} = sprintf('  %-40s  p = %.2e', 'Event rate (KS test)',      p5);
L{end+1} = sprintf('  %-40s  p = %.2e', 'Mean dF/F (KS test)',       p6);

for k = 1:numel(L)
    fprintf('%s\n', L{k});
    fprintf(fid, '%s\n', L{k});
end
fclose(fid);

elapsed = toc;
fprintf('\nAll 7 figures + summary saved to: %s\n', out_dir);
fprintf('Total elapsed: %.1f s (%.1f min)\n', elapsed, elapsed/60);

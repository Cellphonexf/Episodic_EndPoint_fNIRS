%% fNIRS data analysis
%% Required to run in MATLAB version R2017b
%% Required to load Homer3 processed data files
%% Programmed by Feng Xiao (updated on 2025.11.6)
clear all,
clc,
%% Parameter settings
samplingrate = 5.85;
Channels = 1:34;
MediCond = [2 3]; %2 for present-focused exercise; 3 for endpoint-focused exercise
Mindfulness_first_subjs = [1 3 7 9 11 13 15 17 19 29 31 33 35 37 39];
Endpoint_first_subjs = [4 6 8 10 12 14 16 18 20 21 22 23 24 25 26 27 28 30 32 34 36 38 40];
subj = [Mindfulness_first_subjs, Endpoint_first_subjs];
sclConc = 1e6; %convert Conc from Molar to uMolar
%% Data analysis
cd fNIRS_preprocessed_data\derivatives\homer\
interval1 = [0 390]; %in second, for calculating means
ses_mindfulness_hbo = [];
ses_endpoint_hbo = [];
for i_subj = subj
load(num2str(i_subj), '-mat')
channelnum=size(output.dc.dataTimeSeries,3)/34;
setlength=channelnum*34; 
t=output.dcAvg.time;
temp_chPrune = cell2mat(output.misc.mlActAuto);
chPrune = temp_chPrune(1:34,3); %1:valid channels; 0: bad channels
subj_mindfulness_hbo = [];
subj_endpoint_hbo = [];
HbO_mindfulness = [];
HbO_endpoint = [];
    for j_Cond = MediCond
        for k_Ch = Channels
            SigHbO = output.dcAvg.dataTimeSeries(:,(j_Cond-1)*setlength*3+(k_Ch-1)*3+1)*sclConc;
            intind = [find(t >= interval1(1),1,'first') find(t <= interval1(end),1,'last')]; %extract data from 20s before stimulus onset to 390s after it
            HbO_bas_mean = mean(SigHbO(1:intind(1)));
            HbO_bas_std = std(SigHbO(1:intind(1)));
            HbO_act = SigHbO(intind(1):intind(end));
            HbO = (HbO_act - HbO_bas_mean) / HbO_bas_std; %calculate the relative activation based on the baseline (z-standardization)
            if chPrune(k_Ch,1) == 1
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, HbO];
                elseif ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_endpoint = [HbO_endpoint, HbO];
                end
            else
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, zeros(size(HbO,1),1)];
                elseif ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_endpoint = [HbO_endpoint, zeros(size(HbO,1),1)];
                end
            end
        end
    end
    subj_mindfulness_hbo = [subj_mindfulness_hbo, HbO_mindfulness];
    subj_endpoint_hbo = [subj_endpoint_hbo, HbO_endpoint];
    
    ses_mindfulness_hbo = [ses_mindfulness_hbo, subj_mindfulness_hbo];
    ses_endpoint_hbo = [ses_endpoint_hbo, subj_endpoint_hbo]; %col: channels * subj
end

clear subj_mindfulness_hbo subj_endpoint_hbo HbO_mindfulness HbO_endpoint

Ses_mindfulness_hbo = [];
Ses_endpoint_hbo = [];
Ses_mindfulness_hbo_se = [];
Ses_endpoint_hbo_se = [];
Ses_channel_hbo_t = ones(size(ses_mindfulness_hbo, 1), 34);
Ses_channel_hbo_p = ones(size(ses_mindfulness_hbo, 1), 34);
for i = Channels
    i_temp_mindfulness_hbo = [];
    i_temp_endpoint_hbo = [];
    i_mean_mindfulness_hbo = [];
    i_mean_endpoint_hbo = [];
    i_se_mindfulness_hbo = [];
    i_se_endpoint_hbo = [];
    ch_temp_t = ones(size(ses_mindfulness_hbo, 1), 1);
    ch_temp_p = ones(size(ses_mindfulness_hbo, 1), 1);
    for j = 1:size(subj, 2)
        temp_mindfulness_hbo = ses_mindfulness_hbo(:, i+size(Channels, 2)*(j-1));
        temp_endpoint_hbo = ses_endpoint_hbo(:, i+size(Channels, 2)*(j-1));
        
        i_temp_mindfulness_hbo = [i_temp_mindfulness_hbo, temp_mindfulness_hbo];
        i_temp_endpoint_hbo = [i_temp_endpoint_hbo, temp_endpoint_hbo];
    end
    i_temp_mindfulness_hbo = i_temp_mindfulness_hbo(:, any(i_temp_mindfulness_hbo)); 
    i_temp_endpoint_hbo = i_temp_endpoint_hbo(:, any(i_temp_endpoint_hbo)); %delete column with all zeros
    
    i_mean_mindfulness_hbo = mean(i_temp_mindfulness_hbo, 2);%average subjects' data for each channel
    i_mean_endpoint_hbo = mean(i_temp_endpoint_hbo, 2); %average subjects' data for each channel
    
    i_se_mindfulness_hbo = std(i_temp_mindfulness_hbo, 0, 2)./sqrt(size(i_temp_mindfulness_hbo, 2)); %calculate standard error
    i_se_endpoint_hbo = std(i_temp_endpoint_hbo, 0, 2)./sqrt(size(i_temp_endpoint_hbo, 2)); %calculate standard error   
    for k = 1:size(ses_mindfulness_hbo,1)
        [h,p_hbo,ci,stats_hbo] = ttest(i_temp_mindfulness_hbo(k,:), i_temp_endpoint_hbo(k,:)); %paired ttest between two exercises for each time point in each channel
        ch_temp_hbo_t(k) = stats_hbo.tstat; %t>0 denote present-focused > endpoint-focused
        ch_temp_hbo_p(k) = p_hbo;
    end
    Ses_mindfulness_hbo = [Ses_mindfulness_hbo, i_mean_mindfulness_hbo];
    Ses_endpoint_hbo = [Ses_endpoint_hbo, i_mean_endpoint_hbo];
    Ses_mindfulness_hbo_se = [Ses_mindfulness_hbo_se, i_se_mindfulness_hbo];
    Ses_endpoint_hbo_se = [Ses_endpoint_hbo_se, i_se_endpoint_hbo];
    Ses_channel_hbo_t(:, i) = ch_temp_hbo_t;
    Ses_channel_hbo_p(:, i) = ch_temp_hbo_p;
end

clear i_temp_mindfulness_hbo
clear i_temp_endpoint_hbo
clear i_mean_mindfulness_hbo 
clear i_mean_endpoint_hbo 
clear i_se_mindfulness_hbo 
clear i_se_endpoint_hbo
clear temp_mindfulness_hbo 
clear temp_endpoint_hbo 
clear ch_temp_hbo_t 
clear ch_temp_hbo_p 
clear h ci p_hbo stats_hbo 

tTest_Ses_mindfulness_hbo = zeros(5,34);
tTest_Ses_endpoint_hbo    = zeros(5,34);
tTest_Ses_contrast_hbo    = zeros(5,34);

pM = nan(1,34); pE = nan(1,34); pC = nan(1,34);   % temporal for uncorrected p
tM = nan(1,34); tE = nan(1,34); tC = nan(1,34);   % temporal for t-values
ciM = nan(2,34); ciE = nan(2,34); ciC = nan(2,34);% temporal for CI

for i = 1:34
    % Present-focused
    temp_col_m = Ses_mindfulness_hbo(:, i);
    [~, p_m, ci_m, stats_m] = ttest(temp_col_m); % one-sample t-test vs 0
    tM(i)      = stats_m.tstat;
    pM(i)      = p_m;
    ciM(:, i)  = ci_m;

    % Endpoint-focused
    temp_col_i = Ses_endpoint_hbo(:, i);
    [~, p_i, ci_i, stats_i] = ttest(temp_col_i); % one-sample t-test vs 0
    tE(i)      = stats_i.tstat;
    pE(i)      = p_i;
    ciE(:, i)  = ci_i;

    % Contrast: Endpoint - Present£¨channel-wise pair ttests£©
    [~, p_c, ci_c, stats_c] = ttest(temp_col_i, temp_col_m);
    tC(i)      = stats_c.tstat;
    pC(i)      = p_c;
    ciC(:, i)  = ci_c;
end

% BH-FDR for multiple comparison corrections (34 channels)
if exist('mafdr','file') == 2
    qM = mafdr(pM, 'BHFDR', true);
    qE = mafdr(pE, 'BHFDR', true);
    qC = mafdr(pC, 'BHFDR', true);
else
    qM = bh_fdr(pM);
    qE = bh_fdr(pE);
    qC = bh_fdr(pC);
end

% Note: line1(t-values); line2(q-values); line3(lower-CI); line4(upper-CI); line5(significance mark) 
alpha_q = 0.001;  % set alpha level of 0.001
tTest_Ses_mindfulness_hbo(1,:) = tM;
tTest_Ses_mindfulness_hbo(2,:) = qM;           % corrected p-values
tTest_Ses_mindfulness_hbo(3:4,:)= ciM;
tTest_Ses_mindfulness_hbo(5,:) = qM < alpha_q; %1 for accept h1; 0 for accept h0

tTest_Ses_endpoint_hbo(1,:) = tE;
tTest_Ses_endpoint_hbo(2,:) = qE;
tTest_Ses_endpoint_hbo(3:4,:)= ciE;
tTest_Ses_endpoint_hbo(5,:) = qE < alpha_q;

tTest_Ses_contrast_hbo(1,:) = tC;
tTest_Ses_contrast_hbo(2,:) = qC;
tTest_Ses_contrast_hbo(3:4,:)= ciC;
tTest_Ses_contrast_hbo(5,:) = qC < alpha_q;
%% Effect size calculation (Cohen's d, time-series level)
hbo_mindfulness_ses = zeros(3,34);
hbo_endpoint_ses    = zeros(3,34);
efs_hbo_contrast    = zeros(1,34);

for i = 1:34
    % ----- Present-focused£¨mindfulness£© -----
    temp_hbo_m = Ses_mindfulness_hbo(:, i);
    m_mean = mean(temp_hbo_m, 'omitnan');
    m_std  = std(temp_hbo_m, 0, 'omitnan');
    m_se   = m_std / sqrt(length(temp_hbo_m));
    hbo_mindfulness_ses(:, i) = [m_mean; m_se; m_std];

    % ----- Endpoint-focused -----
    temp_hbo_e = Ses_endpoint_hbo(:, i);
    e_mean = mean(temp_hbo_e, 'omitnan');
    e_std  = std(temp_hbo_e, 0, 'omitnan');
    e_se   = e_std / sqrt(length(temp_hbo_e));
    hbo_endpoint_ses(:, i) = [e_mean; e_se; e_std];

    % ----- Cohen's d for contrast (time-series level) -----
    % Endpoint - Present across time points
    diff_series = temp_hbo_e - temp_hbo_m;
    pooled_std = sqrt((m_std^2 + e_std^2) / 2);
    efs_hbo_contrast(1, i) = mean(diff_series, 'omitnan') / pooled_std;
end

efs_hbo_endpoint     = hbo_endpoint_ses(1,:) ./ hbo_endpoint_ses(3,:);
efs_hbo_mindfulness  = hbo_mindfulness_ses(1,:) ./ hbo_mindfulness_ses(3,:);

ES_time_series = table((1:34)', ...
    efs_hbo_mindfulness', efs_hbo_endpoint', efs_hbo_contrast', ...
    'VariableNames', {'Channel','d_mindfulness','d_endpoint','d_contrast'});
%% Visualization for overall activation contrast
cd ..\..\..
cd fNIRS_pics\
%% Line plots (M ¡À SE) for selected channels
% Color setting
col_endpoint = [0.80, 0.10, 0.10];  % red£¨Endpoint-focused£©
col_present  = [0.10, 0.10, 0.80];  % blue£¨Present-focused£©
gray_fill    = [0.85, 0.85, 0.85];  % light gray (mental time travel period)

% X-axis: time (s)
if exist('samplingrate','var') && ~isempty(samplingrate)
    Fs = samplingrate;
else
    Fs = 5.85;
end
nT   = size(Ses_mindfulness_hbo, 1);
Tsec = (0:nT-1)./Fs;
x_min = 0; x_max = 400;
x_mask = (Tsec >= x_min) & (Tsec <= x_max);
Tx = Tsec(x_mask);

% Eight channels with significant contrast differences
chs_to_plot = [1 2 7 11 12 16 18 22];
BAs_to_plot = {'BA 11L','BA 11R','BA 46L','BA 10L','BA 10L','BA 9L','BA 45L','BA 21L'};

% Size of pic
figTS = figure('Units','inches','Position',[0 0 7.0 4.0]); 

% SE
drawShaded = @(x, m, se, c) patch([x, fliplr(x)], [ (m-se)', fliplr((m+se)') ], ...
                                  c, 'FaceAlpha', 0.15, 'EdgeColor','none');

% Subplot letter labels
letters = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)'};                              
                              
for k = 1:numel(chs_to_plot)
    ch = chs_to_plot(k);
    ba = BAs_to_plot{k};

    m_pres = Ses_mindfulness_hbo(x_mask, ch);
    se_pres = Ses_mindfulness_hbo_se(x_mask, ch);
    m_end  = Ses_endpoint_hbo(x_mask, ch);
    se_end = Ses_endpoint_hbo_se(x_mask, ch);

    subplot(2,4,k); hold on;
    
    fill([60 390 390 60], [ -20 -20 20 20 ], gray_fill, ...
         'FaceAlpha', 0.25, 'EdgeColor','none');

    drawShaded(Tx, m_end,  se_end,  col_endpoint);
    drawShaded(Tx, m_pres, se_pres, col_present);
    plot(Tx, m_end,  'Color', col_endpoint, 'LineWidth', 0.8);
    plot(Tx, m_pres, 'Color', col_present,  'LineWidth', 0.8);
    plot([x_min x_max], [0 0], '-', 'Color', [0 0 0 0.15], 'LineWidth', 0.4);

    set(gca, 'FontName','Times New Roman','FontSize',7, ...
        'XColor','k','YColor','k','TickLength',[0.01 0.01], 'Box','off');
    xlim([x_min x_max]);
    set(gca, 'XTick', 0:100:400);
    ylim([-20 20]);
    set(gca, 'YTick', -20:10:20);

    title(sprintf('%s  CH%d - %s', letters{k}, ch, ba), 'FontName','Times New Roman', ...
        'Color','k','FontSize',7,'FontWeight','bold');
    if k > 4
        xlabel('Time (s)','FontName','Times New Roman','FontSize',7,'Color','k');
    end
    ylabel('\DeltaHbO (\muM)','FontName','Times New Roman','FontSize',7,'Color','k'); 

    hold off;
end

% Output for pics
print(figTS, 'fNIRS_timeseries_CHs', '-dpdf', '-r600');
%% Resting-state channel connectivity analysis (full time-series)
num_channels = size(Ses_endpoint_hbo, 2);
T_end = size(Ses_endpoint_hbo, 1);
T_pre = size(Ses_mindfulness_hbo, 1);

% Correlation matrix
r_matrix_endpoint    = corr(Ses_endpoint_hbo,    'rows','pairwise');
r_matrix_mindfulness = corr(Ses_mindfulness_hbo, 'rows','pairwise');

r_matrix_endpoint(    r_matrix_endpoint    >=  1) =  1 - eps;
r_matrix_endpoint(    r_matrix_endpoint    <= -1) = -1 + eps;
r_matrix_mindfulness( r_matrix_mindfulness >=  1) =  1 - eps;
r_matrix_mindfulness( r_matrix_mindfulness <= -1) = -1 + eps;

% Fisher-z transformation
z_matrix_endpoint    = 0.5 * log((1 + r_matrix_endpoint)    ./ (1 - r_matrix_endpoint));
z_matrix_mindfulness = 0.5 * log((1 + r_matrix_mindfulness) ./ (1 - r_matrix_mindfulness));

z_matrix_endpoint(   logical(eye(num_channels))) = NaN;
z_matrix_mindfulness(logical(eye(num_channels))) = NaN;

% Descriptive results for global Fisher-z
mean_z_endpoint    = nanmean(z_matrix_endpoint(:));    % 0.54
mean_z_mindfulness = nanmean(z_matrix_mindfulness(:)); % 0.37

% Between-exercise comparisons (global differences)
z_diff = z_matrix_endpoint - z_matrix_mindfulness;           % Ch x Ch
z_diff(logical(eye(num_channels))) = NaN;
mean_diff   = nanmean(z_diff(:));
std_diff    = nanstd(z_diff(:));
n_edges     = num_channels * (num_channels - 1) / 2;
se_diff     = std_diff / sqrt(n_edges);
t_stat      = mean_diff / se_diff;         % 7.51
df          = n_edges - 1;                 % 560
p_value     = 2 * (1 - tcdf(abs(t_stat), df)); % <.001
cohens_d    = mean_diff / std_diff;        % 0.32

% Edge-wise differences with BH-FDR and dual-criteria flag
maskTri = triu(true(num_channels), 1);
idx     = find(maskTri);
[Ch_i, Ch_j] = ind2sub([num_channels num_channels], idx);

z1 = z_matrix_endpoint(idx);      % Endpoint condition (vectorized edges)
z2 = z_matrix_mindfulness(idx);   % Present  condition (vectorized edges)

% Edge-wise z-test under Fisher-z variance
SE_edge     = sqrt( 2 / (nT - 3) );         % scalar
z_stat_edge = (z1 - z2) ./ SE_edge;         % standard normal under H0
p_edge      = 2 * (1 - normcdf(abs(z_stat_edge), 0, 1));  % two-tailed p-values

% BH-FDR 
[p_sorted, ord] = sort(p_edge(:), 'ascend');  % column vector
m     = numel(p_sorted);
ranks = (1:m)';

q_sorted      = p_sorted .* (m ./ ranks);     % raw BH
q_sorted_rev  = cummin( flipud(q_sorted) );   % monotonicity enforcement
q_sorted_adj  = flipud(q_sorted_rev);
q_sorted_adj  = min(q_sorted_adj, 1);

q = zeros(m,1);
q(ord) = q_sorted_adj;
p_FDR_edge = q;   % FDR-adjusted p-values for edges (vector)

% Effect size per edge: Cohen¡¯s d in Fisher-z space (signed)
d_edge = (z1 - z2);
dual_pass = (p_FDR_edge < 1e-3) & (abs(d_edge) > 0.5);

Edges_results = table( ...
    Ch_i(:), Ch_j(:), ...
    z1(:), z2(:), ...
    (z1(:) - z2(:)), z_stat_edge(:), ...
    p_edge(:), p_FDR_edge(:), ...
    abs(d_edge(:)) > 0.5, dual_pass(:), ...
    'VariableNames', { ...
        'Ch_i','Ch_j', ...
        'Z_endpoint','Z_present', ...
        'Z_diff','Z_stat', ...
        'p_raw','p_FDR', ...
        'd_gt_0p5','DualCriteria_pass'});
    
% Degrees per condition (endpoint / present)
invalid_channels = [10, 13, 21, 23, 25, 27, 30, 31, 33, 34];
if ~exist('num_channels','var'); num_channels = size(z_matrix_endpoint,1); end
if exist('invalid_channels','var')
    valid_channels = setdiff(1:num_channels, invalid_channels);
else
    valid_channels = 1:num_channels;
end

% ---- (A) Degree by |Fisher-z| threshold (default 0.5) ----
thr = 0.5;

Zep = z_matrix_endpoint(valid_channels,    valid_channels);
Zpr = z_matrix_mindfulness(valid_channels, valid_channels);

Aep = double(abs(Zep) >= thr);
Apr = double(abs(Zpr) >= thr);

% remove diagonal; keep simple-graph (undirected, no self-loops)
Aep(1:size(Aep,1)+1:end) = 0;
Apr(1:size(Apr,1)+1:end) = 0;

% use upper triangle to avoid double counting, then symmetrize for degree
Aep_u = triu(Aep,1);
Apr_u = triu(Apr,1);

deg_ep_thr = sum(Aep_u + Aep_u.', 2);           % per-channel degree (endpoint)
deg_pr_thr = sum(Apr_u + Apr_u.', 2);           % per-channel degree (present)
E_ep_thr   = sum(deg_ep_thr)/2;                 % total #edges (endpoint)
E_pr_thr   = sum(deg_pr_thr)/2;                 % total #edges (present)

Degree_byThreshold = table( ...
    valid_channels(:), deg_ep_thr(:), deg_pr_thr(:), ...
    'VariableNames', {'Channel','Degree_ep_thr','Degree_pr_thr'});

% ---- (B) Degree by dual-criteria significant edges (FDR + |d|>0.5) ----
dual_mat = false(num_channels);                           % full graph size
for i = 1:height(Edges_results)
    if Edges_results.DualCriteria_pass(i)
        ii = Edges_results.Ch_i(i);
        jj = Edges_results.Ch_j(i);
        if ~(isnan(ii) || isnan(jj))
            dual_mat(ii,jj) = true;
            dual_mat(jj,ii) = true;
        end
    end
end

dual_mat_valid = dual_mat(valid_channels, valid_channels);
dual_mat_valid(1:size(dual_mat_valid,1)+1:end) = 0;

% per-channel degree on dual-criteria significant graph
deg_dual = sum(dual_mat_valid, 2);
E_dual   = sum(deg_dual)/2;                              % total #edges (dual-criteria)

Degree_byDual = table( ...
    valid_channels(:), deg_dual(:), ...
    'VariableNames', {'Channel','Degree_dualSig'});

% ---- quick text summary in console ----
fprintf('[Threshold |z|>=%.2f] Endpoint edges: %d, Present edges: %d\n', thr, round(E_ep_thr), round(E_pr_thr));
fprintf('[Dual-criteria significant] Edges: %d\n', round(E_dual));  
%% Functional Connectivity (FC) visualization
invalid_channels = [10, 13, 21, 23, 25, 27, 30, 31, 33, 34];
threshold = 0.5; % cut-off for visualizing strong FC
valid_channels = setdiff(1:num_channels, invalid_channels);

% ---------- 1. Endpoint-focused FC ----------
z_ep_tri = tril(z_matrix_endpoint, -1);
z_ep_tri(isnan(z_ep_tri)) = NaN;
z_ep_valid = z_ep_tri(valid_channels, valid_channels);

FC_ep = figure;
set(FC_ep, 'Units','inches','Position',[0,0,4,4]);
imagesc(z_ep_valid);
colormap(autumn);
colorbar;
caxis([-1 1]);
axis square;

ax = gca;
ax.XColor='k'; ax.YColor='k';
ax.TickLength=[0.01 0.01];
xticks(1:numel(valid_channels));
yticks(1:numel(valid_channels));
xticklabels(valid_channels);
yticklabels(valid_channels);
xlabel('Channels','FontName','Times New Roman','FontSize',7,'Color','k');
ylabel('Channels','FontName','Times New Roman','FontSize',7,'Color','k');
title('Functional connectivity: Endpoint-focused','FontName','Times New Roman',...
      'FontSize',7,'FontWeight','bold','Color','k');
set(gca,'FontName','Times New Roman','FontSize',7,'Box','off');
print(FC_ep,'FC_endpoint','-dpdf','-r600');

% Save .edge file (thresholded)
z_ep_sym = z_ep_valid + z_ep_valid';
z_ep_sym(abs(z_ep_sym)<threshold)=0;
dlmwrite('FC_endpoint.edge', z_ep_sym, 'delimiter','\t');

% ---------- 2. Present-focused FC ----------
z_pr_tri = tril(z_matrix_mindfulness, -1);
z_pr_tri(isnan(z_pr_tri)) = NaN;
z_pr_valid = z_pr_tri(valid_channels, valid_channels);

FC_pr = figure;
set(FC_pr, 'Units','inches','Position',[0,0,4,4]);
imagesc(z_pr_valid);
colormap(autumn);
colorbar;
caxis([-1 1]);
axis square;

ax = gca;
ax.XColor='k'; ax.YColor='k';
ax.TickLength=[0.01 0.01];
xticks(1:numel(valid_channels));
yticks(1:numel(valid_channels));
xticklabels(valid_channels);
yticklabels(valid_channels);
xlabel('Channels','FontName','Times New Roman','FontSize',7,'Color','k');
ylabel('Channels','FontName','Times New Roman','FontSize',7,'Color','k');
title('Functional connectivity: Present-focused','FontName','Times New Roman',...
      'FontSize',7,'FontWeight','bold','Color','k');
set(gca,'FontName','Times New Roman','FontSize',7,'Box','off');
print(FC_pr,'FC_present','-dpdf','-r600');

z_pr_sym = z_pr_valid + z_pr_valid';
z_pr_sym(abs(z_pr_sym)<threshold)=0;
dlmwrite('FC_present.edge', z_pr_sym, 'delimiter','\t');

% ---------- 3. Overlay significant edges (dual-criteria passed) ----------
dual_mat_sig = false(num_channels);
color_mat = zeros(num_channels); % 0 = no diff, +1 = endpoint stronger, -1 = present stronger

for i = 1:height(Edges_results)
    if Edges_results.DualCriteria_pass(i)
        ci = Edges_results.Ch_i(i);
        cj = Edges_results.Ch_j(i);
        diff_val = Edges_results.Z_diff(i);
        dual_mat_sig(ci, cj) = true;
        dual_mat_sig(cj, ci) = true;
        % mark direction
        if diff_val > 0
            color_mat(ci, cj) = 1;  % endpoint stronger
            color_mat(cj, ci) = 1;
        elseif diff_val < 0
            color_mat(ci, cj) = -1; % present stronger
            color_mat(cj, ci) = -1;
        end
    end
end

% restrict to valid channels
color_valid = color_mat(valid_channels, valid_channels);

% visualization
FC_sig = figure;
set(FC_sig,'Units','inches','Position',[0,0,4,4]);

imagesc(color_valid);
colormap([0 0 1; 1 1 1; 1 0 0]); % blue-white-red: blue=present>endpoint, red=endpoint>present
caxis([-1 1]);
axis square;

xticks(1:numel(valid_channels));
yticks(1:numel(valid_channels));
xticklabels(valid_channels);
yticklabels(valid_channels);
xlabel('Channels','FontName','Times New Roman','FontSize',7,'Color','k');
ylabel('Channels','FontName','Times New Roman','FontSize',7,'Color','k');
title('Significant edges between two exercises', ...
      'FontName','Times New Roman','FontSize',7,'FontWeight','bold','Color','k');
set(gca,'FontName','Times New Roman','FontSize',7,'Box','off');

print(FC_sig,'FC_significant_edges','-dpdf','-r600');
%% MNI node file generation
[num_data, txt_data, raw_data] = xlsread('FC table.xlsx', 'mni');
mni_data = raw_data;
mni_coords = cell2mat(mni_data(:, 1:3));
channel_labels = cell2mat(mni_data(:, 4));
color = ones(18, 1); %1 for prefrontal lobe's channels
color(19:24) = 2; %2 for temporal lobe's channels
degree_endpoint = sum(abs(z_ep_valid) > 0, 2); %calculate the FC counts for each channel
degree_present = sum(abs(z_pr_valid) > 0, 2); %calculate the FC counts for each channel

endpoint_node = [mni_coords, color, degree_endpoint, channel_labels];
present_node = [mni_coords, color, degree_present, channel_labels];

filename_node = 'nodes_endpoint.node';
dlmwrite(filename_node, endpoint_node, 'delimiter', '\t');
filename_node = 'nodes_present.node';
dlmwrite(filename_node, present_node, 'delimiter', '\t');
%% Hemodynamic change extraction for each participant
%% For the correlation analysis with delay discounting change
%% Endpoint-focused
[num_timepoints, num_total_columns] = size(ses_endpoint_hbo);
num_subjects = num_total_columns / num_channels;
mean_activation_endpoint = zeros(num_subjects, num_channels); %subj*channel
for subj = 1:num_subjects
    start_col = (subj - 1) * num_channels + 1;
    end_col = subj * num_channels;         
    subj_data = ses_endpoint_hbo(:, start_col:end_col); 
    mean_activation_endpoint(subj, :) = mean(subj_data, 1); % mean extraction
end
%% Present-focused
[num_timepoints, num_total_columns] = size(ses_mindfulness_hbo);
num_subjects = num_total_columns / num_channels;
mean_activation_mindfulness = zeros(num_subjects, num_channels); %subj*channel
for subj = 1:num_subjects
    start_col = (subj - 1) * num_channels + 1;
    end_col = subj * num_channels;         
    subj_data = ses_mindfulness_hbo(:, start_col:end_col); 
    mean_activation_mindfulness(subj, :) = mean(subj_data, 1); % mean extraction
end
%% Export mean activation results to Excel
channel_labels = arrayfun(@(x) sprintf('CH%d', x), 1:num_channels, 'UniformOutput', false);
subject_labels = arrayfun(@(x) sprintf('S%02d', x), 1:size(mean_activation_endpoint,1), 'UniformOutput', false);

T_endpoint = array2table(mean_activation_endpoint, 'VariableNames', channel_labels, 'RowNames', subject_labels);
T_mindfulness = array2table(mean_activation_mindfulness, 'VariableNames', channel_labels, 'RowNames', subject_labels);

output_file = 'Mean_HbO_Activation.xlsx';
writetable(T_endpoint, output_file, 'Sheet', 'Endpoint_focused', 'WriteRowNames', true);
writetable(T_mindfulness, output_file, 'Sheet', 'Present_focused', 'WriteRowNames', true);
fprintf('Done! Mean activation data exported successfully to "%s"\n', output_file);
%% fNIRS data analysis (for participants without meditation experiences)
%% Required to run in MATLAB version R2017b
%% Required to load Homer3 processed data files
%% Programmed by Feng Xiao (updated on 2025.11.7)
clear all,
clc,
%% Parameter settings
samplingrate = 5.85;
Channels = 1:34;
MediCond = [2 3]; %2 for present-focused exercise; 3 for endpoint-focused exercise
Mindfulness_first_subjs = [1 3 7 9 11 13 17 19 31 33 35 37 39];
Endpoint_first_subjs = [4 6 8 10 14 16 20 23 26 30 32 34 36 38 40];
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
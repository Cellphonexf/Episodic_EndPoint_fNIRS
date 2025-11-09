%% fNIRS data visualization
%% Programmed by Feng Xiao (updated on 2025.11.6)
clear all,
clc,
%% Load fNIRS_visualization.xlsx
[~, sheetNames] = xlsfinfo('fNIRS_visualization.xlsx');

for i = 1:numel(sheetNames)
    [data, ~, ~] = xlsread('fNIRS_visualization.xlsx', sheetNames{i});
    assignin('base', sheetNames{i}, data);
end
df = 37;
%% Exercise contrast
mni = mni_front;
HbO2 = contrast_front;
Hb = ref_front; %because EasyTopo is required to input both the HbO2 and Hb for visualization, so we used the sham-data (i.e., all data points were -1) as reference
save('input_contrast_data.mat', 'mni', 'HbO2', 'Hb', 'df');
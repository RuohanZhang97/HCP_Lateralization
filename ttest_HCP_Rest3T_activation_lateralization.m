clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
load(fullfile(datapath,'data/HCP_WM_avg_rest_3T_signal.mat'));
load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
savefolder = fullfile(datapath,'results_v3_new');
if ~exist(savefolder,'dir')
    mkdir(savefolder);
end
tbl_savefolder = fullfile(datapath,'results_v3_new','tbl');
if ~exist(tbl_savefolder,'dir')
    mkdir(tbl_savefolder);
end
tblname1 = fullfile(tbl_savefolder,'ttest_HCP_WM_rest_3T_signal_lateralization.xlsx');
tblname2 = fullfile(tbl_savefolder,'ttest_HCP_WM_rest_3T_signal_lateralization_signif_Bon.xlsx');
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis/code'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
rest_cov = [cov.Sex cov.Age cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_Rest3T];

%% ttest between conditions
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
BrainRegionID = cell2mat(LabelID(1:NumRegion,1));
Bon_threshold = 0.05/NumRegion;

[lia,locb] = ismember(cov.SubID,rest_3T_SubID);
rest_3T_signal = rest_3T_signal(:,locb);
NumSub = length(rest_3T_signal);

ttest_results.task = 'Rest';
pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
for m = 1:NumRegion
    tmp_mean_signal = [rest_3T_signal(m,:) rest_3T_signal(m+180,:)]';
    tmp_cov = [rest_cov;rest_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x2,x1); % right - left
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x2,x1,'paired');
end
ttest_results.tVal = tVal;
ttest_results.Cohen_d = Cohen_d;
ttest_results.pVal = pVal;
ttest_results.NumSignifRegion = sum(ttest_results.pVal < Bon_threshold);
ttest_results.NumSub = NumSub;

T = table;
T.BrainRegion = BrainRegionLabel;
T.BrainRegionID = BrainRegionID;
T.tVal = ttest_results.tVal;
T.Cohen_d = ttest_results.Cohen_d;
T.pVal = ttest_results.pVal;
writetable(T,tblname1);
T = sortrows(T,'pVal');
ttest_results.signif_results = T(T.pVal < Bon_threshold,:);
writetable(ttest_results.signif_results,tblname2);

savefile = fullfile(savefolder,'ttest_HCP_Rest3T_activation_lateralization.mat');
save(savefile,'ttest_results');


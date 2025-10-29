clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
files = dir(fullfile(datapath,'data','sorted_selected_HCP_task_block_avg_mean_signal_HCPMMP360'));
files(1:2) = [];
load(fullfile(files(2).folder,files(2).name),'f0bk_faces_mean_signal','f0bk_faces_mean_signal_SubID','f0bk_body_mean_signal','f0bk_body_mean_signal_SubID', ...
    'f0bk_places_mean_signal','f0bk_places_mean_signal_SubID','f0bk_tools_mean_signal','f0bk_tools_mean_signal_SubID',...
    'f0bk_faces_mean_baseline','f0bk_faces_mean_baseline_SubID','f0bk_body_mean_baseline','f0bk_body_mean_baseline_SubID', ...
    'f0bk_places_mean_baseline','f0bk_places_mean_baseline_SubID','f0bk_tools_mean_baseline','f0bk_tools_mean_baseline_SubID');
f0bk_faces = f0bk_faces_mean_signal - f0bk_faces_mean_baseline;
f0bk_places = f0bk_places_mean_signal - f0bk_places_mean_baseline;
f0bk_body = f0bk_body_mean_signal - f0bk_body_mean_baseline;
f0bk_tools = f0bk_tools_mean_signal - f0bk_tools_mean_baseline;

load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
savefolder = fullfile(datapath,'results_v3_new');
if ~exist(savefolder,'dir')
    mkdir(savefolder);
end
tbl_savefolder = fullfile(datapath,'results_v3_new','tbl');
if ~exist(tbl_savefolder,'dir')
    mkdir(tbl_savefolder);
end
tblname1 = fullfile(tbl_savefolder,'ttest_HCP_WM_block_signal_lateralization_between_conditions.xlsx');
tblname2 = fullfile(tbl_savefolder,'ttest_HCP_WM_block_signal_lateralization_between_conditions_signif_Bon.xlsx');
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis/code'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
task_cov = [cov.Age cov.Sex cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_WM];

%% ttest on R minus L between conditions
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
BrainRegionID = cell2mat(LabelID(1:NumRegion,1));
Bon_threshold = 0.05/NumRegion;

ttest_results.task = 'WM';
ttest_results.stimuli = 'f0bk_places_vs_f0bk_faces';
tmp_f0bk_activation_SubID = f0bk_faces_mean_signal_SubID;
[lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, tmp_f0bk_activation_SubID, 'UniformOutput', false)));

tmp_R_f0bk_places = f0bk_places(181:360,locb);
tmp_L_f0bk_places = f0bk_places(1:180,locb);
tmp_RL_f0bk_places = tmp_R_f0bk_places - tmp_L_f0bk_places;

tmp_R_f0bk_faces = f0bk_faces(181:360,locb);
tmp_L_f0bk_faces = f0bk_faces(1:180,locb);
tmp_RL_f0bk_faces = tmp_R_f0bk_faces - tmp_L_f0bk_faces;

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_RL_f0bk_places(m,:) tmp_RL_f0bk_faces(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - faces
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
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
T.places_R = mean(tmp_R_f0bk_places,2);
T.places_L = mean(tmp_L_f0bk_places,2);
T.faces_R = mean(tmp_R_f0bk_faces,2);
T.faces_L = mean(tmp_L_f0bk_faces,2);
writetable(T,tblname1,'Sheet',ttest_results.stimuli);
T = sortrows(T,'pVal');
ttest_results.signif_results = T(T.pVal < Bon_threshold,:);
writetable(ttest_results.signif_results,tblname2,'Sheet',ttest_results.stimuli);

savefile = fullfile(savefolder,'ttest_HCP_WM_block_activation_lateralization_between_conditions.mat');
save(savefile,'ttest_results');


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
tblname1 = fullfile(tbl_savefolder,'ttest_HCP_WM_block_signal_lateralization.xlsx');
tblname2 = fullfile(tbl_savefolder,'ttest_HCP_WM_block_signal_lateralization_signif_Bon.xlsx');
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis/code'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
task_cov = [cov.Age cov.Sex cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_WM];

%% ttest between conditions
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
BrainRegionID = cell2mat(LabelID(1:NumRegion,1));
Bon_threshold = 0.05/NumRegion;
stimuli = {'f0bk_body','f0bk_faces','f0bk_places','f0bk_tools'};

for i = 1:length(stimuli)
    disp(i)
    ttest_results(i).task = 'WM';
    ttest_results(i).stimuli = stimuli{i};
    tmp_f0bk_activation_SubID = eval([stimuli{i},'_mean_signal_SubID']);
    [lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, tmp_f0bk_activation_SubID, 'UniformOutput', false)));
    tmp_f0bk_activation = eval(stimuli{i});
    tmp_f0bk_activation = tmp_f0bk_activation(:,locb);

    pVal = zeros(NumRegion,1);
    tVal = zeros(NumRegion,1);
    NumSub = size(tmp_f0bk_activation,2);
    for m = 1:NumRegion
        tmp_mean_signal = [tmp_f0bk_activation(m,:) tmp_f0bk_activation(m+180,:)]';
        tmp_cov = [task_cov;task_cov];
        [~,~,r] = regress(tmp_mean_signal,tmp_cov);
        x1 = r(1:NumSub,1);
        x2 = r(NumSub+1:end,1);
        [~,p,~,stats] = ttest(x2,x1); % right - left
        pVal(m,1) = p;
        tVal(m,1) = stats.tstat;
        Cohen_d(m,1) = computeCohen_d(x2,x1,'paired');
    end
    ttest_results(i).tVal = tVal;
    ttest_results(i).Cohen_d = Cohen_d;
    ttest_results(i).pVal = pVal;
    ttest_results(i).NumSignifRegion = sum(ttest_results(i).pVal < Bon_threshold);
    ttest_results(i).NumSub = NumSub;

    T = table;
    T.BrainRegion = BrainRegionLabel;
    T.BrainRegionID = BrainRegionID;
    T.tVal = ttest_results(i).tVal;
    T.Cohen_d = ttest_results(i).Cohen_d;
    T.pVal = ttest_results(i).pVal;
    writetable(T,tblname1,'Sheet',ttest_results(i).stimuli);
    T = sortrows(T,'pVal');
    ttest_results(i).signif_results = T(T.pVal < Bon_threshold,:);
    writetable(ttest_results(i).signif_results,tblname2,'Sheet',ttest_results(i).stimuli);
end

savefile = fullfile(savefolder,'ttest_HCP_WM_block_activation_lateralization.mat');
save(savefile,'ttest_results');


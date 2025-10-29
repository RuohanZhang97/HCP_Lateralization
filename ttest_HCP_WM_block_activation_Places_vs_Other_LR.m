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
% tbl_savefolder = fullfile(datapath,'results_v3_new','tbl');
% if ~exist(tbl_savefolder,'dir')
%     mkdir(tbl_savefolder);
% end
% tblname1 = fullfile(tbl_savefolder,'ttest_HCP_WM_block_signal_lateralization_places_minus_faces.xlsx');
% tblname2 = fullfile(tbl_savefolder,'ttest_HCP_WM_block_signal_lateralization_places_minus_faces_signif_Bon.xlsx');
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis/code'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
task_cov = [cov.Age cov.Sex cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_WM];

%% ttest on Places vs Other stimuli in each hemisphere
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
BrainRegionID = cell2mat(LabelID(1:NumRegion,1));
Bon_threshold = 0.05/NumRegion;

tmp_f0bk_activation_SubID = f0bk_places_mean_signal_SubID;
[lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, tmp_f0bk_activation_SubID, 'UniformOutput', false)));

% places vs faces Left
ttest_results(1).task = 'WM';
ttest_results(1).stimuli = 'places_vs_faces_L';

tmp_L_f0bk_places = f0bk_places(1:180,locb);
tmp_L_f0bk_faces = f0bk_faces(1:180,locb);

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_L_f0bk_places(m,:) tmp_L_f0bk_faces(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - others
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
end

ttest_results(1).tVal = tVal;
ttest_results(1).Cohen_d = Cohen_d;
ttest_results(1).pVal = pVal;
ttest_results(1).NumSignifRegion = sum(ttest_results(1).pVal < Bon_threshold);
ttest_results(1).NumSub = NumSub;

% T = table;
% T.BrainRegion = BrainRegionLabel;
% T.BrainRegionID = BrainRegionID;
% T.tVal = ttest_results.tVal;
% T.Cohen_d = ttest_results.Cohen_d;
% T.pVal = ttest_results.pVal;
% T.places_L = mean(tmp_L_f0bk_places,2);
% T.faces_L = mean(tmp_L_f0bk_faces,2);
% writetable(T,tblname1,'Sheet',ttest_results.stimuli);
% T = sortrows(T,'pVal');
% ttest_results.signif_results = T(T.pVal < Bon_threshold,:);
% writetable(ttest_results.signif_results,tblname2,'Sheet',ttest_results.stimuli);

% places vs faces Right
ttest_results(2).task = 'WM';
ttest_results(2).stimuli = 'places_vs_faces_R';

tmp_R_f0bk_places = f0bk_places(181:360,locb);
tmp_R_f0bk_faces = f0bk_faces(181:360,locb);

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_R_f0bk_places(m,:) tmp_R_f0bk_faces(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - others
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
end

ttest_results(2).tVal = tVal;
ttest_results(2).Cohen_d = Cohen_d;
ttest_results(2).pVal = pVal;
ttest_results(2).NumSignifRegion = sum(ttest_results(2).pVal < Bon_threshold);
ttest_results(2).NumSub = NumSub;


% places vs body Left
ttest_results(3).task = 'WM';
ttest_results(3).stimuli = 'places_vs_body_L';

tmp_L_f0bk_places = f0bk_places(1:180,locb);
tmp_L_f0bk_body = f0bk_body(1:180,locb);

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_L_f0bk_places(m,:) tmp_L_f0bk_body(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - others
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
end

ttest_results(3).tVal = tVal;
ttest_results(3).Cohen_d = Cohen_d;
ttest_results(3).pVal = pVal;
ttest_results(3).NumSignifRegion = sum(ttest_results(3).pVal < Bon_threshold);
ttest_results(3).NumSub = NumSub;

% places vs body Right
ttest_results(4).task = 'WM';
ttest_results(4).stimuli = 'places_vs_body_R';

tmp_R_f0bk_places = f0bk_places(181:360,locb);
tmp_R_f0bk_body = f0bk_body(181:360,locb);

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_R_f0bk_places(m,:) tmp_R_f0bk_body(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - others
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
end

ttest_results(4).tVal = tVal;
ttest_results(4).Cohen_d = Cohen_d;
ttest_results(4).pVal = pVal;
ttest_results(4).NumSignifRegion = sum(ttest_results(4).pVal < Bon_threshold);
ttest_results(4).NumSub = NumSub;


% places vs tools Left
ttest_results(5).task = 'WM';
ttest_results(5).stimuli = 'places_vs_tools_L';

tmp_L_f0bk_places = f0bk_places(1:180,locb);
tmp_L_f0bk_tools = f0bk_tools(1:180,locb);

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_L_f0bk_places(m,:) tmp_L_f0bk_tools(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - others
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
end

ttest_results(5).tVal = tVal;
ttest_results(5).Cohen_d = Cohen_d;
ttest_results(5).pVal = pVal;
ttest_results(5).NumSignifRegion = sum(ttest_results(5).pVal < Bon_threshold);
ttest_results(5).NumSub = NumSub;

% places vs tools Right
ttest_results(6).task = 'WM';
ttest_results(6).stimuli = 'places_vs_tools_R';

tmp_R_f0bk_places = f0bk_places(181:360,locb);
tmp_R_f0bk_tools = f0bk_tools(181:360,locb);

pVal = zeros(NumRegion,1);
tVal = zeros(NumRegion,1);
NumSub = length(locb);
for m = 1:NumRegion
    tmp_mean_signal = [tmp_R_f0bk_places(m,:) tmp_R_f0bk_tools(m,:)]';
    tmp_cov = [task_cov;task_cov];
    [~,~,r] = regress(tmp_mean_signal,tmp_cov);
    x1 = r(1:NumSub,1);
    x2 = r(NumSub+1:end,1);
    [~,p,~,stats] = ttest(x1,x2); % places - others
    pVal(m,1) = p;
    tVal(m,1) = stats.tstat;
    Cohen_d(m,1) = computeCohen_d(x1,x2,'paired');
end

ttest_results(6).tVal = tVal;
ttest_results(6).Cohen_d = Cohen_d;
ttest_results(6).pVal = pVal;
ttest_results(6).NumSignifRegion = sum(ttest_results(6).pVal < Bon_threshold);
ttest_results(6).NumSub = NumSub;


savefile = fullfile(savefolder,'ttest_HCP_WM_block_activation_places_vs_others_LR.mat');
save(savefile,'ttest_results');

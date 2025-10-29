clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
files = dir(fullfile(datapath,'data','sorted_selected_HCP_task_block_rest_avg_FC_HCPMMP360_V2'));
files(1:2) = [];
task_files = files;
load(fullfile(task_files(1).folder,task_files(1).name),'f0bk_body_FC','f0bk_body_FC_SubID','f0bk_faces_FC','f0bk_faces_FC_SubID',...
    'f0bk_places_FC','f0bk_places_FC_SubID','f0bk_tools_FC','f0bk_tools_FC_SubID');
load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
savefolder = fullfile(datapath,'results_v3_new');
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/code/UKB_FC_analysis'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
basic_cov = [cov.Age cov.Sex cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_WM];

%% t-Test on FCs for places vs others
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
NumLink = NumRegion*(NumRegion-1)/2;
Bon_threshold = 0.05/NumLink;
task_names = {'WM','WM','WM','WM'};
stimuli = {'f0bk_body','f0bk_faces','f0bk_tools'};

[lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, f0bk_places_FC_SubID, 'UniformOutput', false)));
tmp_f0bk_places_FC = f0bk_places_FC(locb,:,:);
tmp_avg_places_FC_mtx = squeeze(mean(tmp_f0bk_places_FC,1));

% places vs others Left
for i = 1:length(stimuli)
    disp(i)

    tmp_f0bk_FC_SubID = eval([stimuli{i},'_FC_SubID']);
    [lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, tmp_f0bk_FC_SubID, 'UniformOutput', false)));
    tmp_f0bk_FC = eval([stimuli{i},'_FC']);
    tmp_f0bk_FC = tmp_f0bk_FC(locb,:,:);

    ttest_results_L(i).task = 'WM';
    ttest_results_L(i).stimuli = ['f0bk_places_vs_' stimuli{i}];

    task_cov = [basic_cov; basic_cov];
    pVal = zeros(NumRegion,NumRegion);
    tVal = zeros(NumRegion,NumRegion);
    Cohen_d = zeros(NumRegion,NumRegion);
    for m = 1:NumRegion-1
        for n = m+1:NumRegion
            tmp_FC = [tmp_f0bk_places_FC(:,m,n); tmp_f0bk_FC(:,m,n)];
            [~,~,r] = regress(tmp_FC,task_cov);
            x1 = r(1:length(tmp_FC)/2,1);
            x2 = r(length(tmp_FC)/2+1:end,1);
            [~,p,~,stats] = ttest(x1,x2); % places - others
            pVal(m,n) = p;
            tVal(m,n) = stats.tstat;
            Cohen_d(m,n) = computeCohen_d(x1,x2,'paired');
        end
    end
    tmp_avg_FC_mtx = squeeze(mean(tmp_f0bk_FC,1));
    
    ttest_results_L(i).tVal = tVal;
    ttest_results_L(i).Cohen_d = Cohen_d;
    ttest_results_L(i).pVal = pVal;
    ttest_results_L(i).FDR = Vector2FCtriu(mafdr(FCtriu2Vector(pVal,1),'BHFDR',true),1);
    ttest_results_L(i).NumSignifLinkBon = sum(FCtriu2Vector(ttest_results_L(i).pVal,1) < Bon_threshold);
    ttest_results_L(i).NumSignifLinkFDR = sum(ttest_results_L(i).FDR < 0.05);
    ttest_results_L(i).NumSub = size(tmp_f0bk_FC,1);
    ttest_results_L(i).FCdiff = tmp_avg_places_FC_mtx(1:180,1:180) - tmp_avg_FC_mtx(1:180,1:180); % places minus others
end

% places vs others Right
for i = 1:length(stimuli)
    disp(i)

    tmp_f0bk_FC_SubID = eval([stimuli{i},'_FC_SubID']);
    [lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, tmp_f0bk_FC_SubID, 'UniformOutput', false)));
    tmp_f0bk_FC = eval([stimuli{i},'_FC']);
    tmp_f0bk_FC = tmp_f0bk_FC(locb,:,:);

    ttest_results_R(i).task = 'WM';
    ttest_results_R(i).stimuli = ['f0bk_places_vs_' stimuli{i}];

    task_cov = [basic_cov; basic_cov];
    pVal = zeros(NumRegion,NumRegion);
    tVal = zeros(NumRegion,NumRegion);
    Cohen_d = zeros(NumRegion,NumRegion);
    for m = 1:NumRegion-1
        for n = m+1:NumRegion
            tmp_FC = [tmp_f0bk_places_FC(:,m+180,n+180); tmp_f0bk_FC(:,m+180,n+180)];
            [~,~,r] = regress(tmp_FC,task_cov);
            x1 = r(1:length(tmp_FC)/2,1);
            x2 = r(length(tmp_FC)/2+1:end,1);
            [~,p,~,stats] = ttest(x1,x2); % places - others
            pVal(m,n) = p;
            tVal(m,n) = stats.tstat;
            Cohen_d(m,n) = computeCohen_d(x1,x2,'paired');
        end
    end
    tmp_avg_FC_mtx = squeeze(mean(tmp_f0bk_FC,1));
    
    ttest_results_R(i).tVal = tVal;
    ttest_results_R(i).Cohen_d = Cohen_d;
    ttest_results_R(i).pVal = pVal;
    ttest_results_R(i).FDR = Vector2FCtriu(mafdr(FCtriu2Vector(pVal,1),'BHFDR',true),1);
    ttest_results_R(i).NumSignifLinkBon = sum(FCtriu2Vector(ttest_results_R(i).pVal,1) < Bon_threshold);
    ttest_results_R(i).NumSignifLinkFDR = sum(ttest_results_R(i).FDR < 0.05);
    ttest_results_R(i).NumSub = size(tmp_f0bk_FC,1);
    ttest_results_R(i).FCdiff = tmp_avg_places_FC_mtx(181:360,181:360) - tmp_avg_FC_mtx(181:360,181:360); % places minus others
end

savefile = fullfile(savefolder,'ttest_HCP_WM_FC_places_vs_others_LR.mat');
save(savefile,'ttest_results_L','ttest_results_R');



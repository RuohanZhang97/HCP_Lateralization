% This code is written by Ruohan Zhang on 14th Feb 2025 for the
% system-relvant FC lateralization tests and then FDR and Bon corrections.
clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
savefolder = fullfile(datapath,'results_v3_new');
load(fullfile(savefolder,'ttest_HCP_WM_FC_lateralization.mat'));

fig_savefolder1 = fullfile(savefolder,'fig','HCP_WM_FC_lateralization_HippSys','Bon');
if ~exist(fig_savefolder1,'dir')
    mkdir(fig_savefolder1);
end
fig_savefolder2 = fullfile(savefolder,'fig','HCP_WM_FC_lateralization_HippSys','FDR');
if ~exist(fig_savefolder2,'dir')
    mkdir(fig_savefolder2);
end
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/code/UKB_FC_analysis'));

load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
task_names = {'WM','WM','WM','WM'};
stimuli = {'f0bk_body','f0bk_faces','f0bk_places','f0bk_tools'};

% The following code is to select System * 180 regions, reduce repeated links, and do FDR and Bon corrections
MasterInd = [80:87 14:17 127 131 128:130 132];
SysBrainLabel = BrainRegionLabel(MasterInd,1);
BrainRegionLabel1 = BrainRegionLabel(1:90,1);
BrainRegionLabel2 = BrainRegionLabel(91:180,1);

for s = 1:length(ttest_results)
    RightMinusLeftFC = ttest_results(s).RL_FCdiff;
    pVal = ttest_results(s).pVal;
    Sys_RightMinusLeftFC = RightMinusLeftFC(MasterInd,:);

    selected_pVal = [];
    index_pairs = []; % record (i,j)
    for i = 1:length(MasterInd)
        for j = 1:NumRegion
            if MasterInd(i) ~= j
                index_pairs = [index_pairs; MasterInd(i),j]; % record index
            end
        end
    end
    sorted_pairs = index_pairs;
    sorted_pairs = sort(index_pairs, 2);
    [~, unique_idx] = unique(sorted_pairs, 'rows');
    unique_pairs = index_pairs(unique_idx, :);

    pVal = pVal' + pVal;
    for i = 1:length(unique_pairs)
        selected_pVal = [selected_pVal; pVal(unique_pairs(i,1), unique_pairs(i,2))];
    end

    % FDR correction
    FDR_selected_pVal = mafdr(selected_pVal, 'BHFDR', true);
    FDR_selected_pVal = FDR_selected_pVal .* (FDR_selected_pVal < 0.05);
    FDR_matrix = [unique_pairs FDR_selected_pVal];
    FDR_matrix(FDR_matrix(:,3) == 0, :) = [];

    % Bonferroni correction
    Bon_threshold = 0.05/length(selected_pVal);
    Bon_selected_pVal = selected_pVal .* (selected_pVal < Bon_threshold);
    Bon_matrix = [unique_pairs Bon_selected_pVal];
    Bon_matrix(Bon_matrix(:,3) == 0, :) = [];

    % FDR correction matrix for Sys*180 brain regions
    idx_matrix = nan(length(MasterInd), NumRegion);
    for i = 1:size(FDR_matrix,1)
        row_id = FDR_matrix(i,1);
        col_id = FDR_matrix(i,2);
        row_index = find(MasterInd == row_id);
        if ~isempty(row_index)
            idx_matrix(row_index, col_id) = 1;
        end
    end
    FDR_Sys_RightMinusLeftFC = idx_matrix .* Sys_RightMinusLeftFC;

    % Bon correction matrix for Sys*180 brain regions
    idx_matrix = nan(length(MasterInd), NumRegion);
    for i = 1:size(Bon_matrix,1)
        row_id = Bon_matrix(i,1);
        col_id = Bon_matrix(i,2);
        row_index = find(MasterInd == row_id);
        if ~isempty(row_index)
            idx_matrix(row_index, col_id) = 1;
        end
    end
    Bon_Sys_RightMinusLeftFC = idx_matrix .* Sys_RightMinusLeftFC;

    % Plot these two matrix FDR_Sys_RightMinusLeftFC and Bon_Sys_RightMinusLeftFC
    % FDR corrections first 90 regions
    FDR_Sys_RightMinusLeftFC1 = FDR_Sys_RightMinusLeftFC(:,1:90);
    figure(1)
    fig1 = imagesc(FDR_Sys_RightMinusLeftFC1);
    set(fig1,'AlphaData',~isnan(FDR_Sys_RightMinusLeftFC1));
    set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

    % FDR corrections last 90 regions
    FDR_Sys_RightMinusLeftFC2 = FDR_Sys_RightMinusLeftFC(:,91:180);
    figure(2)
    fig2 = imagesc(FDR_Sys_RightMinusLeftFC2);
    set(fig2,'AlphaData',~isnan(FDR_Sys_RightMinusLeftFC2));
    set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

    % Bon corrections first 90 regions
    Bon_Sys_RightMinusLeftFC1 = Bon_Sys_RightMinusLeftFC(:,1:90);
    figure(3)
    fig1 = imagesc(Bon_Sys_RightMinusLeftFC1);
    set(fig1,'AlphaData',~isnan(Bon_Sys_RightMinusLeftFC1));
    set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

    % Bon corrections last 90 regions
    Bon_Sys_RightMinusLeftFC2 = Bon_Sys_RightMinusLeftFC(:,91:180);
    figure(4)
    fig2 = imagesc(Bon_Sys_RightMinusLeftFC2);
    set(fig2,'AlphaData',~isnan(Bon_Sys_RightMinusLeftFC2));
    set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

end


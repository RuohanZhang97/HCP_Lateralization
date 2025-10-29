% This code is written by Ruohan Zhang on 13 Sep 2025 for the
% system-relvant FC lateralization tests and then FDR and Bon corrections.
clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
savefolder = fullfile(datapath,'results_v3_new');
load(fullfile(savefolder,'ttest_HCP_WM_FC_places_vs_others_LR.mat'));

fig_savefolder1 = fullfile(savefolder,'fig','HCP_WM_FC_places_vs_others_Visual','Bon');
if ~exist(fig_savefolder1,'dir')
    mkdir(fig_savefolder1);
end
fig_savefolder2 = fullfile(savefolder,'fig','HCP_WM_FC_places_vs_others_Visual','FDR');
if ~exist(fig_savefolder2,'dir')
    mkdir(fig_savefolder2);
end
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/code/UKB_FC_analysis'));

load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
NumRegion = 180;
BrainRegionLabel = LabelID(1:NumRegion,4);

% The following code is to select System * 180 regions, reduce repeated links, and do FDR and Bon corrections
MasterInd = [1:26 40 46 47 62:66 84:95 106:113]; % visual  % [80:87 14:17 127 131 128:130 132]; % hipp
SysBrainLabel = BrainRegionLabel(MasterInd,1);
BrainRegionLabel1 = BrainRegionLabel(1:90,1);
BrainRegionLabel2 = BrainRegionLabel(91:180,1);

%% FC differences for places vs others Left hemisphere
for s = 1:length(ttest_results_L)
    disp(s)
    FCdiff_L = ttest_results_L(s).FCdiff;
    pVal = ttest_results_L(s).pVal;
    Sys_FCdiff_L = FCdiff_L(MasterInd,:);

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

%     % Bonferroni correction
%     Bon_threshold = 0.05/length(selected_pVal);
%     Bon_selected_pVal = selected_pVal .* (selected_pVal < Bon_threshold);
%     Bon_matrix = [unique_pairs Bon_selected_pVal];
%     Bon_matrix(Bon_matrix(:,3) == 0, :) = [];

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
    FDR_Sys_FCdiff_L = idx_matrix .* Sys_FCdiff_L;

%     % Bon correction matrix for Sys*180 brain regions
%     idx_matrix = nan(length(MasterInd), NumRegion);
%     for i = 1:size(Bon_matrix,1)
%         row_id = Bon_matrix(i,1);
%         col_id = Bon_matrix(i,2);
%         row_index = find(MasterInd == row_id);
%         if ~isempty(row_index)
%             idx_matrix(row_index, col_id) = 1;
%         end
%     end
%     Bon_Sys_FC_diff_L = idx_matrix .* Sys_FCdiff_L;

    % Plot FDR_Sys_FCdiff
    % FDR corrections first 90 regions
    FDR_Sys_FCdiff1_L = FDR_Sys_FCdiff_L(:,1:90);
    figure(1)
    fig1 = imagesc(FDR_Sys_FCdiff1_L);
    set(fig1,'AlphaData',~isnan(FDR_Sys_FCdiff1_L));
    set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

    % FDR corrections last 90 regions
    FDR_Sys_FCdiff2_L = FDR_Sys_FCdiff_L(:,91:180);
    figure(2)
    fig2 = imagesc(FDR_Sys_FCdiff2_L);
    set(fig2,'AlphaData',~isnan(FDR_Sys_FCdiff2_L));
    set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

%     % Bon corrections first 90 regions
%     Bon_Sys_FCdiff1_L = Bon_Sys_FC_diff_L(:,1:90);
%     figure(3)
%     fig1 = imagesc(Bon_Sys_FCdiff1_L);
%     set(fig1,'AlphaData',~isnan(Bon_Sys_FCdiff1_L));
%     set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1),'FontSize',10);
%     xtickangle(90)
%     set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
%     set(gca,'TickLength',[0 0]);
%     Xtick_pos = 1.5:NumRegion;
%     Ytick_pos = 1.5:length(MasterInd);
%     xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
%     yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);
% 
%     % Bon corrections last 90 regions
%     Bon_Sys_FCdiff2_L = Bon_Sys_FC_diff_L(:,91:180);
%     figure(4)
%     fig2 = imagesc(Bon_Sys_FCdiff2_L);
%     set(fig2,'AlphaData',~isnan(Bon_Sys_FCdiff2_L));
%     set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2),'FontSize',10);
%     xtickangle(90)
%     set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
%     set(gca,'TickLength',[0 0]);
%     Xtick_pos = 1.5:NumRegion;
%     Ytick_pos = 1.5:length(MasterInd);
%     xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
%     yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

end

%% FC differences for places vs others Right hemisphere
for s = 1:length(ttest_results_R)
    disp(s)
    FCdiff_R = ttest_results_R(s).FCdiff;
    pVal = ttest_results_R(s).pVal;
    Sys_FCdiff_R = FCdiff_R(MasterInd,:);

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

%     % Bonferroni correction
%     Bon_threshold = 0.05/length(selected_pVal);
%     Bon_selected_pVal = selected_pVal .* (selected_pVal < Bon_threshold);
%     Bon_matrix = [unique_pairs Bon_selected_pVal];
%     Bon_matrix(Bon_matrix(:,3) == 0, :) = [];

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
    FDR_Sys_FCdiff_R = idx_matrix .* Sys_FCdiff_R;

%     % Bon correction matrix for Sys*180 brain regions
%     idx_matrix = nan(length(MasterInd), NumRegion);
%     for i = 1:size(Bon_matrix,1)
%         row_id = Bon_matrix(i,1);
%         col_id = Bon_matrix(i,2);
%         row_index = find(MasterInd == row_id);
%         if ~isempty(row_index)
%             idx_matrix(row_index, col_id) = 1;
%         end
%     end
%     Bon_Sys_FC_diff_R = idx_matrix .* Sys_FCdiff_R;

    % Plot FDR_Sys_FCdiff
    % FDR corrections first 90 regions
    FDR_Sys_FCdiff1_R = FDR_Sys_FCdiff_R(:,1:90);
    figure(1)
    fig1 = imagesc(FDR_Sys_FCdiff1_R);
    set(fig1,'AlphaData',~isnan(FDR_Sys_FCdiff1_R));
    set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);
    
    % FDR corrections last 90 regions
    FDR_Sys_FCdiff2_R = FDR_Sys_FCdiff_R(:,91:180);
    figure(2)
    fig2 = imagesc(FDR_Sys_FCdiff2_R);
    set(fig2,'AlphaData',~isnan(FDR_Sys_FCdiff2_R));
    set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2),'FontSize',10);
    xtickangle(90)
    set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
    set(gca,'TickLength',[0 0]);
    Xtick_pos = 1.5:NumRegion;
    Ytick_pos = 1.5:length(MasterInd);
    xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
    yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

%     % Bon corrections first 90 regions
%     Bon_Sys_FCdiff1_R = Bon_Sys_FC_diff_R(:,1:90);
%     figure(3)
%     fig1 = imagesc(Bon_Sys_FCdiff1_R);
%     set(fig1,'AlphaData',~isnan(Bon_Sys_FCdiff1_R));
%     set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1),'FontSize',10);
%     xtickangle(90)
%     set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
%     set(gca,'TickLength',[0 0]);
%     Xtick_pos = 1.5:NumRegion;
%     Ytick_pos = 1.5:length(MasterInd);
%     xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
%     yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);
% 
%     % Bon corrections last 90 regions
%     Bon_Sys_FCdiff2_R = Bon_Sys_FC_diff_R(:,91:180);
%     figure(4)
%     fig2 = imagesc(Bon_Sys_FCdiff2_R);
%     set(fig2,'AlphaData',~isnan(Bon_Sys_FCdiff2_R));
%     set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2),'FontSize',10);
%     xtickangle(90)
%     set(gca, 'YTickLabel',SysBrainLabel,'YTick',1:length(SysBrainLabel),'FontSize',10);
%     set(gca,'TickLength',[0 0]);
%     Xtick_pos = 1.5:NumRegion;
%     Ytick_pos = 1.5:length(MasterInd);
%     xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
%     yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);

end


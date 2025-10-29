clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
files = dir(fullfile(datapath,'data','sorted_selected_HCP_task_block_rest_avg_FC_HCP_MMP360'));
files(1:2) = [];
rest_file = files(6);
load(fullfile(rest_file(1).folder,rest_file(1).name));
load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
savefolder = fullfile(datapath,'results_v3_new');
fig_savefolder1 = fullfile(savefolder,'fig','HCP_Rest3T_FC_lateralization','Bon');
if ~exist(fig_savefolder1,'dir')
    mkdir(fig_savefolder1);
end
fig_savefolder2 = fullfile(savefolder,'fig','HCP_Rest3T_FC_lateralization','FDR');
if ~exist(fig_savefolder2,'dir')
    mkdir(fig_savefolder2);
end
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/code/UKB_FC_analysis'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
basic_cov = [cov.Age cov.Sex cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_Rest3T];

%% t-Test on lateralization in FCs
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
NumLink = NumRegion*(NumRegion-1)/2;
Bon_threshold = 0.05/NumLink;

[lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, FC_SubID, 'UniformOutput', false)));
FC = FC(locb,:,:);

final_cov = [basic_cov; basic_cov];
pVal = zeros(NumRegion,NumRegion);
tVal = zeros(NumRegion,NumRegion);
Cohen_d = zeros(NumRegion,NumRegion);
for m = 1:NumRegion-1
    for n = m+1:NumRegion
        tmp_FC = [FC(:,m,n); FC(:,m+180,n+180)];
        [~,~,r] = regress(tmp_FC,final_cov);
        x1 = r(1:length(tmp_FC)/2,1);
        x2 = r(length(tmp_FC)/2+1:end,1);
        [~,p,~,stats] = ttest(x2,x1); % right - left
        pVal(m,n) = p;
        tVal(m,n) = stats.tstat;
        Cohen_d(m,n) = computeCohen_d(x2,x1,'paired');
    end
end
tmp_avg_FC_mtx = squeeze(mean(FC,1));
ttest_results.stimuli = 'Rest3T';
ttest_results.tVal = tVal;
ttest_results.Cohen_d = Cohen_d;
ttest_results.pVal = pVal;
ttest_results.FDR = Vector2FCtriu(mafdr(FCtriu2Vector(pVal,1),'BHFDR',true),1);
ttest_results.NumSignifLinkBon = sum(FCtriu2Vector(ttest_results.pVal,1) < Bon_threshold);
ttest_results.NumSignifLinkFDR = sum(ttest_results.FDR < 0.05);
ttest_results.NumSub = size(FC,1);
ttest_results.RL_FCdiff = tmp_avg_FC_mtx(181:360,181:360)-tmp_avg_FC_mtx(1:180,1:180); % Right minus Left

savefile = fullfile(savefolder,'ttest_HCP_Rest3T_FC_lateralization.mat');
save(savefile,'ttest_results');

%% Figures for FC maps (Bon and FDR corrections)
fig_scale_value = 1.2;
Cohen_d = ttest_results.Cohen_d;

%% Bonferroni correction
upper_d_matrix = triu(Cohen_d .* (ttest_results.pVal < Bon_threshold),1);
lower_d_matrix = Cohen_d';
[top_dval,ind_top_dval] = sort(FCtriu2Vector(abs(upper_d_matrix),1),'descend');
PropThreshold = 0.4; % 0.3; % 0.2; % 1.0;
ind_top_Bon_signif = zeros(NumLink,1);
ind_top_Bon_signif(ind_top_dval(1:round(sum(top_dval > 0)*PropThreshold))) = 1;
ind_top_signif_matrix = Vector2FCtriu(ind_top_Bon_signif,1); % transform the indices for the significant links into 0 or 1 (1 represents significant links)
d_matrix = upper_d_matrix .* ind_top_signif_matrix + lower_d_matrix;
d_matrix(d_matrix == 0) = NaN;

close all;
figure(1)
fig = imagesc(d_matrix);
set(fig,'AlphaData',~isnan(d_matrix));
set(gca,'XTickLabel',BrainRegionLabel,'XTick',1:length(BrainRegionLabel),'FontSize',2.5);
xtickangle(90)
set(gca, 'YTickLabel',BrainRegionLabel,'YTick',1:length(BrainRegionLabel),'FontSize',2.5);
set(gca,'TickLength',[0 0]);
Xtick_pos = 1.5:NumRegion;
Ytick_pos = 1.5:NumRegion;
xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);
FigName = fullfile(fig_savefolder1,strcat('HCP_WM_',ttest_results.stimuli,'_FC_lateralization_prop',num2str(PropThreshold),'.png'));

[rows, column] = size(d_matrix);
% indicate brain system with black lines
Region_Index=[1 4 10 17 26 31 40 47 52 59 67 79 87 95 100 110 120 133 149 158 167];
xline(Region_Index+0.5,'-k','LineWidth',0.9);
yline(Region_Index+0.5,'-k','LineWidth',0.9);

cb = colorbar;
cb.Label.String = 'Cohen''s d';
cb.Label.FontSize = 9;
cb.FontSize = 8;

colormap('jet');
caxis([-fig_scale_value fig_scale_value]);
cb.Ticks = [-fig_scale_value:0.4:fig_scale_value];

f = gcf;
set(f,'PaperUnits','centimeters','PaperPosition',[0 0 18 18]);
print(figure(1),FigName,'-dpng','-r600');

%% FDR correction
upper_d_matrix = triu(Cohen_d .* (ttest_results.FDR < 0.05),1);
lower_d_matrix = Cohen_d';
[top_dval,ind_top_dval] = sort(FCtriu2Vector(abs(upper_d_matrix),1),'descend');
PropThreshold = 0.4; % 0.3; % 0.2; % 1.0;
ind_top_Bon_signif = zeros(NumLink,1);
ind_top_Bon_signif(ind_top_dval(1:round(sum(top_dval > 0)*PropThreshold))) = 1;
ind_top_signif_matrix = Vector2FCtriu(ind_top_Bon_signif,1); % transform the indices for the significant links into 0 or 1 (1 represents significant links)
d_matrix = upper_d_matrix .* ind_top_signif_matrix + lower_d_matrix;
d_matrix(d_matrix == 0) = NaN;

close all;
figure(1)
fig = imagesc(d_matrix);
set(fig,'AlphaData',~isnan(d_matrix));
set(gca,'XTickLabel',BrainRegionLabel,'XTick',1:length(BrainRegionLabel),'FontSize',2.5);
xtickangle(90)
set(gca, 'YTickLabel',BrainRegionLabel,'YTick',1:length(BrainRegionLabel),'FontSize',2.5);
set(gca,'TickLength',[0 0]);
Xtick_pos = 1.5:NumRegion;
Ytick_pos = 1.5:NumRegion;
xline(Xtick_pos,'-','Color',[0.5 0.5 0.5]);
yline(Ytick_pos,'-','Color',[0.5 0.5 0.5]);
FigName = fullfile(fig_savefolder2,strcat('HCP_WM_',ttest_results.stimuli,'_FC_lateralization_prop',num2str(PropThreshold),'.png'));

[rows, column] = size(d_matrix);
% indicate brain system with black lines
Region_Index=[1 4 10 17 26 31 40 47 52 59 67 79 87 95 100 110 120 133 149 158 167];
xline(Region_Index+0.5,'-k','LineWidth',0.9);
yline(Region_Index+0.5,'-k','LineWidth',0.9);

cb = colorbar;
cb.Label.String = 'Cohen''s d';
cb.Label.FontSize = 9;
cb.FontSize = 8;

colormap('jet');
caxis([-fig_scale_value fig_scale_value]);
cb.Ticks = [-fig_scale_value:0.4:fig_scale_value];

f = gcf;
set(f,'PaperUnits','centimeters','PaperPosition',[0 0 18 18]);
print(figure(1),FigName,'-dpng','-r600');


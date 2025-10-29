%% ========================================================================
%  HCP WM Activation Lateralization (with color-coded significance markers)
%  Right > Left = red, Left > Right = blue
% ========================================================================

clear; clc; close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis/';
load(fullfile(datapath,'data/sorted_selected_HCP_task_block_avg_mean_signal_rawBOLD_HCPMMP360/HCP_WM_block_mean_signal_rawBOLD_HCPMMP360.mat'))

% match subjects by handedness
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
[~,locb] = ismember(cov.SubID,avg_sub_block_signal_SubID);
f0bk_faces = f0bk_faces_mean_signal(:,locb);
f0bk_places = f0bk_places_mean_signal(:,locb);
f0bk_body = f0bk_body_mean_signal(:,locb);
f0bk_tools = f0bk_tools_mean_signal(:,locb);

% load region labels
load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat');
BrainRegionLabel1 = strrep(LabelID(1:90,4),' L','');
BrainRegionLabel2 = strrep(LabelID(91:180,4),' L','');

savefolder = fullfile(datapath,'results_v3_new');
fig_savefolder = fullfile(savefolder,'fig','HCP_WM_activation_3T_rest_LH_RH');
if ~exist(fig_savefolder,'dir'); mkdir(fig_savefolder); end

load(fullfile(savefolder,'ttest_HCP_WM_block_activation_lateralization.mat'));

width = 1800; % figure width
height = 600; % figure height
minVal = 0.4*10^4; % y axis min
maxVal = 1.8*10^4; % y axis max

%% ======================== f0bk_faces ========================
faces_LH = mean(f0bk_faces(1:180,:),2);
faces_RH = mean(f0bk_faces(181:360,:),2);

SignifBonBrainID = sort(ttest_results(2).signif_results.BrainRegionID);
first90 = SignifBonBrainID(SignifBonBrainID <= 90);
last90  = SignifBonBrainID(SignifBonBrainID > 90)-90;

% -------- First 90 regions --------
close all; figure(1)
plot(faces_LH(1:90),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(faces_RH(1:90),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91]);
ylabel('BOLD Signal','FontSize',20);
legend({'Left','Right'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:90);

xtickangle(90);
x_star=first90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,faces_LH,faces_RH)

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_faces_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% -------- Last 90 regions --------
close all; figure(1)
plot(faces_LH(91:180),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(faces_RH(91:180),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167]-90;
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91]);
ylabel('BOLD Signal','FontSize',20);
legend({'Left','Right'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:90);

xtickangle(90);
x_star=last90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,faces_LH(91:180),faces_RH(91:180))

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_faces_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

%% ======================== f0bk_places ========================
places_LH = mean(f0bk_places(1:180,:),2);
places_RH = mean(f0bk_places(181:360,:),2);

SignifBonBrainID = sort(ttest_results(3).signif_results.BrainRegionID);
first90 = SignifBonBrainID(SignifBonBrainID <= 90);
last90  = SignifBonBrainID(SignifBonBrainID > 90)-90;

% First 90
close all; figure(1)
plot(places_LH(1:90),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(places_RH(1:90),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91]);
ylabel('BOLD Signal','FontSize',20)
legend({'Left','Right'},'FontSize',14)
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:90);

xtickangle(90);
x_star=first90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,places_LH,places_RH)

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_places_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% Last 90
close all; figure(1)
plot(places_LH(91:180),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(places_RH(91:180),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167]-90;
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91])
ylabel('BOLD Signal','FontSize',20)
legend({'Left','Right'},'FontSize',14)
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:90);

xtickangle(90);
x_star=last90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,places_LH(91:180),places_RH(91:180))

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_places_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

%% ======================== f0bk_body ========================
body_LH = mean(f0bk_body(1:180,:),2);
body_RH = mean(f0bk_body(181:360,:),2);

SignifBonBrainID = sort(ttest_results(1).signif_results.BrainRegionID);
first90 = SignifBonBrainID(SignifBonBrainID <= 90);
last90  = SignifBonBrainID(SignifBonBrainID > 90)-90;

% First 90
close all; figure(1)
plot(body_LH(1:90),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(body_RH(1:90),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91])
ylabel('BOLD Signal','FontSize',20)
legend({'Left','Right'},'FontSize',14)
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:90);

xtickangle(90);
x_star=first90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,body_LH,body_RH)

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_body_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);


% Last 90
close all; figure(1)
plot(body_LH(91:180),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(body_RH(91:180),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167]-90;
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91])
ylabel('BOLD Signal','FontSize',20)
legend({'Left','Right'},'FontSize',14)
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:90); 

xtickangle(90);
x_star=last90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,body_LH(91:180),body_RH(91:180))

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_body_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

%% ======================== f0bk_tools ========================
tools_LH = mean(f0bk_tools(1:180,:),2);
tools_RH = mean(f0bk_tools(181:360,:),2);

SignifBonBrainID = sort(ttest_results(4).signif_results.BrainRegionID);
first90 = SignifBonBrainID(SignifBonBrainID <= 90);
last90  = SignifBonBrainID(SignifBonBrainID > 90)-90;

% First 90
close all; figure(1)
plot(tools_LH(1:90),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(tools_RH(1:90),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91])
ylabel('BOLD Signal','FontSize',20)
legend({'Left','Right'},'FontSize',14)
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:90); 

xtickangle(90);
x_star=first90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,tools_LH,tools_RH)

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_tools_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% Last 90
close all; figure(1)
plot(tools_LH(91:180),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14); hold on
plot(tools_RH(91:180),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167]-90;
xline(Region_Index+0.5,'-k','LineWidth',0.5)
ylim([minVal maxVal]); xlim([0 91])
ylabel('BOLD Signal','FontSize',20)
legend({'Left','Right'},'FontSize',14)
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:90);

xtickangle(90);
x_star=last90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,tools_LH(91:180),tools_RH(91:180))

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_tools_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

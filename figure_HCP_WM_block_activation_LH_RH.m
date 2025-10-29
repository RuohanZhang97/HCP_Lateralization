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

load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
[lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, f0bk_faces_mean_signal_SubID, 'UniformOutput', false)));
f0bk_faces = f0bk_faces(:,locb);
f0bk_places = f0bk_places(:,locb);
f0bk_body = f0bk_body(:,locb);
f0bk_tools = f0bk_tools(:,locb);

load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
BrainRegionLabel1 = strrep(LabelID(1:90,4),' L','');
BrainRegionLabel2 = strrep(LabelID(91:180,4),' L','');

savefolder = fullfile(datapath,'results_v3_new');
fig_savefolder = fullfile(savefolder,'fig','HCP_WM_activation_3T_rest_LH_RH');
if ~exist(fig_savefolder,'dir')
    mkdir(fig_savefolder);
end

load(fullfile(savefolder,'ttest_HCP_WM_block_activation_lateralization.mat'));

width = 1800; % figure width
height = 600; % figure height
minVal = -60; % y axis min
maxVal = 160; % y axis max

%% f0bk_faces
faces_LH = mean(f0bk_faces(1:180,:),2);
faces_RH = mean(f0bk_faces(181:360,:),2);

SignifBonBrainID = sort(ttest_results(2).signif_results.BrainRegionID);
first90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID <= 90);
last90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID > 90)-90;

% first 90 brain regions
close all
figure(1)
plot(faces_LH(1:90,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(faces_RH(1:90,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1));
xtickangle(90);
% add * symbol to mark significant brain regions after Bon correction
x_star = first90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_faces_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% last 90 brain regions
close all
figure(1)
plot(faces_LH(91:180,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(faces_RH(91:180,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167] - 90;
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2));
xtickangle(90);
% add * symbol to mark significant brain regions after Bon correction
x_star = last90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_faces_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

%% f0bk_places
places_LH = mean(f0bk_places(1:180,:),2);
places_RH = mean(f0bk_places(181:360,:),2);

SignifBonBrainID = sort(ttest_results(3).signif_results.BrainRegionID);
first90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID <= 90);
last90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID > 90)-90;

% first 90 brain regions
close all
figure(1)
plot(places_LH(1:90,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(places_RH(1:90,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1));
xtickangle(90);

% add * symbol to mark significant brain regions after Bon correction
x_star = first90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_places_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% last 90 brain regions
close all
figure(1)
plot(places_LH(91:180,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(places_RH(91:180,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167] - 90;
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2));
xtickangle(90);
% add * symbol to mark significant brain regions after Bon correction
x_star = last90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_places_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

%% f0bk_body
body_LH = mean(f0bk_body(1:180,:),2);
body_RH = mean(f0bk_body(181:360,:),2);

SignifBonBrainID = sort(ttest_results(1).signif_results.BrainRegionID);
first90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID <= 90);
last90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID > 90)-90;

% first 90 brain regions
close all
figure(1)
plot(body_LH(1:90,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(body_RH(1:90,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1));
xtickangle(90);
% add * symbol to mark significant brain regions after Bon correction
x_star = first90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_body_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% last 90 brain regions
close all
figure(1)
plot(body_LH(91:180,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(body_RH(91:180,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167] - 90;
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2));
xtickangle(90);
% add * symbol to mark significant brain regions after Bon correction
x_star = last90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_body_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

%% f0bk_tools
tools_LH = mean(f0bk_tools(1:180,:),2);
tools_RH = mean(f0bk_tools(181:360,:),2);

SignifBonBrainID = sort(ttest_results(4).signif_results.BrainRegionID);
first90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID <= 90);
last90_SignifBonBrainID = SignifBonBrainID(SignifBonBrainID > 90)-90;

% first 90 brain regions
close all
figure(1)
plot(tools_LH(1:90,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(tools_RH(1:90,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1));
xtickangle(90);
box off
% add * symbol to mark significant brain regions after Bon correction
x_star = first90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_tools_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);

% last 90 brain regions
close all
figure(1)
plot(tools_LH(91:180,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(tools_RH(91:180,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167] - 90;
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('Activation','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2));
xtickangle(90);
box off
% add * symbol to mark significant brain regions after Bon correction
x_star = last90_SignifBonBrainID;
y_star = repmat(-50,length(x_star),1); % the min value for the plot is -60
text(x_star, y_star, '*', 'FontSize', 18, 'Color', 'red', 'HorizontalAlignment', 'center');

ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_f0bk_tools_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);


clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
load(fullfile(datapath,'data/HCP_WM_avg_rest_3T_signal.mat'));
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
[lia,locb] = ismember(cov.SubID,rest_3T_SubID);
rest_3T_signal = rest_3T_signal(:,locb);

load('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/data/HCPex_LabelID.mat'); % HCPex labels for brain regions
BrainRegionLabel1 = strrep(LabelID(1:90,4),' L','');
BrainRegionLabel2 = strrep(LabelID(91:180,4),' L','');

savefolder = fullfile(datapath,'results_v3_new');
fig_savefolder = fullfile(savefolder,'fig','HCP_WM_activation_3T_rest_LH_RH');
if ~exist(fig_savefolder,'dir')
    mkdir(fig_savefolder);
end

load(fullfile(savefolder,'ttest_HCP_Rest3T_activation_lateralization.mat'));

width = 1800; % figure width
height = 600; % figure height
minVal = 0.4*10^4; % y axis min
maxVal = 1.8*10^4; % y axis max

%% 3T resting state
rest_LH = mean(rest_3T_signal(1:180,:),2);
rest_RH = mean(rest_3T_signal(181:360,:),2);

SignifBonBrainID = sort(ttest_results.signif_results.BrainRegionID);
first90 = SignifBonBrainID(SignifBonBrainID <= 90);
last90 = SignifBonBrainID(SignifBonBrainID > 90)-90;

% first 90 brain regions
close all
figure(1)
plot(rest_LH(1:90,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(rest_RH(1:90,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[1 4 10 17 26 31 40 47 52 57 66 79 87];
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('BOLD Signal','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel1,'XTick',1:length(BrainRegionLabel1));
xtickangle(90);
x_star=first90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,rest_LH,rest_RH)


box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_rest3T_first90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);


% last 90 brain regions
close all
figure(1)
plot(rest_LH(91:180,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'Marker','.','MarkerSize',14)
hold on
plot(rest_RH(91:180,1),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'Marker','.','MarkerSize',14)
Region_Index=[95 100 110 120 133 149 158 167] - 90;
xline(Region_Index+0.5,'-k','LineWidth',0.5);
ylim([minVal maxVal])
xlim([0 91])
ylabel('BOLD Signal','FontSize',20);
legend({'Left Hemisphere', 'Right Hemisphere'},'FontSize',14);
set(gca,'XTickLabel',BrainRegionLabel2,'XTick',1:length(BrainRegionLabel2));
xtickangle(90);
x_star=first90; y_star = repmat(0.5*10^4,length(x_star),1); % the min value for the plot is -60
plotColoredStars(x_star,y_star,rest_LH,rest_RH)

box off
ax = gca; % 获取当前轴对象
ax.TickLength = [0.005 0.005]; % 设置主刻度和次刻度线的长度
ax.XAxis.FontSize = 12; % 修改 x 轴字体大小
ax.YAxis.FontSize = 18; % 修改 y 轴字体大小
rectangle('Position', [ax.XLim(1), ax.YLim(1), diff(ax.XLim), diff(ax.YLim)],'EdgeColor', 'k', 'LineWidth', 0.9) % 设置边框颜色和宽度
set(gcf,'Position',[1 1 width height]);
FigName = fullfile(fig_savefolder,'HCP_WM_rest3T_last90regions.png');
exportgraphics(gca, FigName, 'Resolution', 600);


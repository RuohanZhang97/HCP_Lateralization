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
fig_savefolder1 = fullfile(savefolder,'fig','HCP_WM_FC_lateralization','Bon');
if ~exist(fig_savefolder1,'dir')
    mkdir(fig_savefolder1);
end
fig_savefolder2 = fullfile(savefolder,'fig','HCP_WM_FC_lateralization','FDR');
if ~exist(fig_savefolder2,'dir')
    mkdir(fig_savefolder2);
end
addpath(genpath('/mnt/disk1/home/ruohan/ScientificProject/UKBdata_analysis/code/UKB_FC_analysis'));

% covariates
load(fullfile(datapath,'data','HCP_cov_HeadMotion_new.mat'));
cov = cov(cov.Handedness >= 40,:);
basic_cov = [cov.Age cov.Sex cov.Drinker_status cov.Smoking_status cov.Education cov.Head_motion_WM];

%% t-Test on lateralization in FCs
NumRegion = 180;
BrainRegionLabel = strrep(LabelID(1:NumRegion,4),' L','');
NumLink = NumRegion*(NumRegion-1)/2;
Bon_threshold = 0.05/NumLink;
task_names = {'WM','WM','WM','WM'};
stimuli = {'f0bk_body','f0bk_faces','f0bk_places','f0bk_tools'};

for i = 1:length(stimuli)
    disp(i)
    ttest_results(i).task = 'WM';
    ttest_results(i).stimuli = stimuli{i};

    tmp_f0bk_FC_SubID = eval([stimuli{i},'_FC_SubID']);
    [lia,locb] = ismember(cov.SubID,cell2mat(cellfun(@str2num, tmp_f0bk_FC_SubID, 'UniformOutput', false)));
    tmp_f0bk_FC = eval([stimuli{i},'_FC']);
    tmp_f0bk_FC = tmp_f0bk_FC(locb,:,:);

    task_cov = [basic_cov; basic_cov];
    pVal = zeros(NumRegion,NumRegion);
    tVal = zeros(NumRegion,NumRegion);
    Cohen_d = zeros(NumRegion,NumRegion);
    for m = 1:NumRegion-1
        for n = m+1:NumRegion
            tmp_FC = [tmp_f0bk_FC(:,m,n); tmp_f0bk_FC(:,m+180,n+180)];
            [~,~,r] = regress(tmp_FC,task_cov);
            x1 = r(1:length(tmp_FC)/2,1);
            x2 = r(length(tmp_FC)/2+1:end,1);
            [~,p,~,stats] = ttest(x2,x1); % right - left
            pVal(m,n) = p;
            tVal(m,n) = stats.tstat;
            Cohen_d(m,n) = computeCohen_d(x2,x1,'paired');
        end
    end
    tmp_avg_FC_mtx = squeeze(mean(eval([stimuli{i} '_FC']),1));
    ttest_results(i).tVal = tVal;
    ttest_results(i).Cohen_d = Cohen_d;
    ttest_results(i).pVal = pVal;
    ttest_results(i).FDR = Vector2FCtriu(mafdr(FCtriu2Vector(pVal,1),'BHFDR',true),1);
    ttest_results(i).NumSignifLinkBon = sum(FCtriu2Vector(ttest_results(i).pVal,1) < Bon_threshold);
    ttest_results(i).NumSignifLinkFDR = sum(ttest_results(i).FDR < 0.05);
    ttest_results(i).NumSub = size(tmp_f0bk_FC,1);
    ttest_results(i).RL_FCdiff = tmp_avg_FC_mtx(181:360,181:360)-tmp_avg_FC_mtx(1:180,1:180); % Right minus Left
end
savefile = fullfile(savefolder,'ttest_HCP_WM_FC_lateralization.mat');
save(savefile,'ttest_results');

%% Figures for FC maps (Bon and FDR corrections) 
% threshold = 1; % Bon correction
threshold = 2; % FDR correction
fig_scale_value = [0.3 0.3 0.5 0.5];
if threshold == 1
    for i = 1:length(ttest_results)
        Cohen_d = ttest_results(i).Cohen_d;
        % Draw d-matrix figure (left hemisphere)
        upper_d_matrix = triu(Cohen_d .* (ttest_results(i).pVal < Bon_threshold),1);
        lower_d_matrix = Cohen_d';
        [top_dval,ind_top_dval] = sort(FCtriu2Vector(abs(upper_d_matrix),1),'descend');
        PropThreshold = 1.0;
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
        FigName = fullfile(fig_savefolder1,strcat('HCP_WM_',ttest_results(i).stimuli,'_FC_lateralization.png'));

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
        caxis([-fig_scale_value(i) fig_scale_value(i)]);
        cb.Ticks = [-fig_scale_value(i):0.1:fig_scale_value(i)];

        f = gcf;
        set(f,'PaperUnits','centimeters','PaperPosition',[0 0 18 18]);
        print(figure(1),FigName,'-dpng','-r600');
    end

else % FDR correction
    for i = 1:length(ttest_results)
        Cohen_d = ttest_results(i).Cohen_d;
        % Draw d-matrix figure (left hemisphere)
        upper_d_matrix = triu(Cohen_d .* (ttest_results(i).FDR < 0.05),1);
        lower_d_matrix = Cohen_d';
        [top_dval,ind_top_dval] = sort(FCtriu2Vector(abs(upper_d_matrix),1),'descend');
        PropThreshold = 1.0;
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
        FigName = fullfile(fig_savefolder2,strcat('HCP_WM_',ttest_results(i).stimuli,'_FC_lateralization.png'));

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
        caxis([-fig_scale_value(i) fig_scale_value(i)]);
        cb.Ticks = [-fig_scale_value(i):0.1:fig_scale_value(i)];

        f = gcf;
        set(f,'PaperUnits','centimeters','PaperPosition',[0 0 18 18]);
        print(figure(1),FigName,'-dpng','-r600');
    end
end


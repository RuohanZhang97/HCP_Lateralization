clear
clc

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis/data/';
load(fullfile(datapath,'selected_HCP_task_block_timeseries_HCP_MMP360/HCP_WM_block_timeseries_HCP_MMP360.mat'));
savefolder = fullfile(datapath,'sorted_selected_HCP_task_block_avg_mean_signal_rawBOLD_HCPMMP360');


for i = 1:length(sub_block_RL)
    disp(i)
    avg_sub_block(i).SubID = sub_block_RL(i).SubID;

    body_points = min(size(sub_block_LR(i).f0bk_body_1,2),size(sub_block_RL(i).f0bk_body_1,2));
    avg_sub_block(i).f0bk_body = (sub_block_LR(i).f0bk_body_1(:,1:body_points) + sub_block_RL(i).f0bk_body_1(:,1:body_points)) ./ 2;

    faces_points = min(size(sub_block_LR(i).f0bk_faces_1,2),size(sub_block_RL(i).f0bk_faces_1,2));
    avg_sub_block(i).f0bk_faces = (sub_block_LR(i).f0bk_faces_1(:,1:faces_points) + sub_block_RL(i).f0bk_faces_1(:,1:faces_points)) ./ 2;
    
    places_points = min(size(sub_block_LR(i).f0bk_places_1,2),size(sub_block_RL(i).f0bk_places_1,2));
    avg_sub_block(i).f0bk_places = (sub_block_LR(i).f0bk_places_1(:,1:places_points) + sub_block_RL(i).f0bk_places_1(:,1:places_points)) ./ 2;

    tools_points = min(size(sub_block_LR(i).f0bk_tools_1,2),size(sub_block_RL(i).f0bk_tools_1,2));
    avg_sub_block(i).f0bk_tools = (sub_block_LR(i).f0bk_tools_1(:,1:tools_points) + sub_block_RL(i).f0bk_tools_1(:,1:tools_points)) ./ 2;

end


for i = 1:length(avg_sub_block)
    disp(i)
    avg_sub_block_signal(i).SubID = avg_sub_block(i).SubID;
    avg_sub_block_signal(i).f0bk_body = mean(avg_sub_block(i).f0bk_body(:,16:end),2); % the first 15 timepoints are prestimuli
    avg_sub_block_signal(i).f0bk_faces = mean(avg_sub_block(i).f0bk_faces(:,16:end),2);
    avg_sub_block_signal(i).f0bk_places = mean(avg_sub_block(i).f0bk_places(:,16:end),2);
    avg_sub_block_signal(i).f0bk_tools = mean(avg_sub_block(i).f0bk_tools(:,16:end),2);

end


avg_sub_block_signal_SubID = str2num(cell2mat({avg_sub_block_signal.SubID}'));
for i = 1:length(avg_sub_block_signal)
    f0bk_body_mean_signal(:,i) = avg_sub_block_signal(i).f0bk_body;
    f0bk_faces_mean_signal(:,i) = avg_sub_block_signal(i).f0bk_faces;
    f0bk_places_mean_signal(:,i) = avg_sub_block_signal(i).f0bk_places;
    f0bk_tools_mean_signal(:,i) = avg_sub_block_signal(i).f0bk_tools;

end

savefile = fullfile(savefolder,'HCP_WM_block_mean_signal_rawBOLD_HCPMMP360.mat');
save(savefile,'avg_sub_block_signal_SubID','f0bk_body_mean_signal','f0bk_faces_mean_signal','f0bk_places_mean_signal','f0bk_tools_mean_signal');



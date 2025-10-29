clear
clc
close all

datapath = '/mnt/disk1/home/ruohan/ScientificProject/HCPdata_analysis';
load(fullfile(datapath,'data/selected_HCP_task_block_timeseries_HCP_MMP360/HCP_WM_block_timeseries_HCP_MMP360.mat'))

sub_block = struct([]);
for i = 1:length(sub_block_LR)
    disp(i)
    sub_block(i).SubID = str2num(sub_block_LR(i).SubID);

    f0bk_body_bin = min(size(sub_block_LR(i).f0bk_body_1,2), size(sub_block_RL(i).f0bk_body_1,2));
    sub_block(i).f0bk_body = (sub_block_LR(i).f0bk_body_1(:,1:f0bk_body_bin) + sub_block_RL(i).f0bk_body_1(:,1:f0bk_body_bin)) / 2;

    f0bk_faces_bin = min(size(sub_block_LR(i).f0bk_faces_1,2), size(sub_block_RL(i).f0bk_faces_1,2));
    sub_block(i).f0bk_faces = (sub_block_LR(i).f0bk_faces_1(:,1:f0bk_faces_bin) + sub_block_RL(i).f0bk_faces_1(:,1:f0bk_faces_bin)) / 2;

    f0bk_places_bin = min(size(sub_block_LR(i).f0bk_places_1,2), size(sub_block_RL(i).f0bk_places_1,2));
    sub_block(i).f0bk_places = (sub_block_LR(i).f0bk_places_1(:,1:f0bk_places_bin) + sub_block_RL(i).f0bk_places_1(:,1:f0bk_places_bin)) / 2;

    f0bk_tools_bin = min(size(sub_block_LR(i).f0bk_tools_1,2), size(sub_block_RL(i).f0bk_tools_1,2));
    sub_block(i).f0bk_tools = (sub_block_LR(i).f0bk_tools_1(:,1:f0bk_tools_bin) + sub_block_RL(i).f0bk_tools_1(:,1:f0bk_tools_bin)) / 2;
end

avg_sub_block = struct([]);
for i = 1:length(sub_block)
    disp(i)
    avg_sub_block(i).SubID = sub_block(i).SubID;
    avg_sub_block(i).f0bk_body = mean(sub_block(i).f0bk_body(:,1:15),"all");
    avg_sub_block(i).f0bk_faces = mean(sub_block(i).f0bk_faces(:,1:15),"all");
    avg_sub_block(i).f0bk_places = mean(sub_block(i).f0bk_places(:,1:15),"all");
    avg_sub_block(i).f0bk_tools = mean(sub_block(i).f0bk_tools(:,1:15),"all");
end


savefile = fullfile(datapath,'data/HCP_WM_block_BOLD_timeseries.mat');
save(savefile,'sub_block','avg_sub_block');





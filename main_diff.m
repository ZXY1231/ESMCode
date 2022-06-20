% clear all;
format long;
% Main function
% Matlab is pass-by-value.
% | Version | Author | Date     | Commit
% | 0.1     | ZhouXY | 20.02.12 | The init version
% | 0.2     | ZhouXY | 20.07.03 | Reconstruct the model for compatbility
% | 0.3     | ZhouXY | 20.07.09 | add gap filling
% | 1.0     | ZhouXY | 20.07.31 | Modify model structure, commit to github
% | 1.1     | ZhouXY | 21.03.03 | For PSM filtering
% | 1.2     | ZhouXY | 21.03.28 | For PSM statistics
% TODO: put the parameters outside function but main function
%test version is used for testing center fitting in time series
%% % Parameters
tic;
frames_path = 'D:\20211017\5nmIgG_5nmIgA_4ms_200fps_1fold_avg10\';

result_path = ["D:\20211017\5nmIgG_5nmIgA_4ms_200fps_1fold_avg10_diff\"
                ];

averageN = 1;
move_step = 1;
%% load images structure and initialize filter operated images, detected particles
imges_struct = dir(append(frames_path ,'*.tif'));
imge_num = length(imges_struct);
shapes = size(imread(append(imges_struct(1).folder, '/' ,imges_struct(1).name)));

%%
tic
leng = floor(imge_num/averageN)*averageN;
for i = 1:leng-2*averageN+1
    image_stack1 = single(imread(append(imges_struct(i).folder, '/' ,imges_struct(i).name)));
    image_stack2 = single(imread(append(imges_struct(i+1).folder, '/' ,imges_struct(i+1).name)));
    stack1_averageN = squeeze(mean(reshape(image_stack1,averageN,[],shapes(1),shapes(2)),1));
    stack2_averageN = squeeze(mean(reshape(image_stack2,averageN,[],shapes(1),shapes(2)),1));
    averageN_diff = single(stack2_averageN - stack1_averageN);
    
    folder_indx = ceil(i/1800000000);
    saveastifffast(averageN_diff, append(result_path(folder_indx), "\",num2str(i,'%06d'), ".tif"));
end
toc






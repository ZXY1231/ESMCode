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
frames_path = 'D:\20211022\532nm_108mW\5nmIgG_5nmIgA_6ms_160fps_1fold_files\';

result_path = ["D:\20211017\5nmIgG_5nmIgA_4ms_200fps_1fold_avg10\" ,...
                ];

averageN = 10;
%%
imges = dir(append(frames_path ,'*.tif'));
imge_num = length(imges);
shapes = size(imread(append(imges(1).folder, '/' ,imges(1).name)));

%%
tic

leng = floor(imge_num/averageN)*averageN;
% leng = 120000;

for i = 1:averageN:leng-2*averageN+1
    image_stack1 = zeros(averageN, shapes(1), shapes(2),'int16');
    
    for j = i:i+averageN-1
        img = imread(append(imges(j).folder, '/' ,imges(j).name));
        image_stack1(mod(j-1,averageN)+1,:,:) = img;
    end
    
    stack1_averageN = squeeze(mean(reshape(image_stack1,averageN,[],shapes(1),shapes(2)),1));

    
    folder_indx = ceil(i/1800000000);
    saveastifffast(single(stack1_averageN), append(result_path(folder_indx), num2str(i,'%06d'), ".tif"));
end
toc






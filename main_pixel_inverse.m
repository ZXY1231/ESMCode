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
frames_path = 'C:\Users\xzhou145\Desktop\20220401\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_060001_110000\';

result_path = ["C:\Users\xzhou145\Desktop\20220401\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_060001_110000_pixel_inverse\" ,...
                ];

%%
all_images = LoadImages(frames_path);% size (#frames,h,w), change images type to single

%%
tic
img_s = size(all_images);
leng = img_s(1);

for i = 1:leng
    image = squeeze(all_images(i,:,:));
    image_inverse = single(-image);
    
    folder_indx = ceil(i/1800000000);
    saveastifffast(image_inverse, append(result_path(folder_indx), num2str(i,'%06d'), ".tif"));
end
toc






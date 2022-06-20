% clear all;
format long;
% Main function

%% % Parameters
tic;
frames_path = 'D:\20211017\5nmIgG_5nmIgA_4ms_200fps_1fold_avg10_diff\';

result_path = ["D:\20211017\5nmIgG_5nmIgA_4ms_200fps_1fold_avg10_diff_img_seletion\"
                ];

averageN = 1;
%% load images structure and initialize filter operated images, detected particles
imges_struct = dir(append(frames_path ,'*.tif'));
imge_num = length(imges_struct);
shapes = size(imread(append(imges_struct(1).folder, '/' ,imges_struct(1).name)));

%%
tic
stdss = zeros(1,imge_num);
for i = 1:imge_num
    image = single(imread(append(imges_struct(i).folder, '/' ,imges_struct(i).name)));
    stdss(i) = std(image,0,'all');
    if std(image,0,'all')>5.2
        continue
    end
 
    folder_indx = ceil(i/1800000000);
    saveastifffast(image, append(result_path(folder_indx), "\",num2str(i,'%06d'), ".tif"));
end
toc

%%
figure(1)
subplot(1,8,1)
plot(stdss)

stdss2 = stdss(stdss<(mean(stdss,'all')+3*std(stdss)));
subplot(1,8,2)
plot(stdss2)

stdss3 = stdss2(stdss2<(mean(stdss2,'all')+3*std(stdss2)));
subplot(1,8,3)
plot(stdss3)

stdss4 = stdss3(stdss3<(mean(stdss3,'all')+3*std(stdss3)));
subplot(1,8,4)
plot(stdss4)

stdss5 = stdss4(stdss4<(mean(stdss4,'all')+3*std(stdss4)));
subplot(1,8,5)
plot(stdss5)

stdss6 = stdss5(stdss5<(mean(stdss5,'all')+3*std(stdss5)));
subplot(1,8,6)
plot(stdss6)

stdss7 = stdss6(stdss6<(mean(stdss6,'all')+3*std(stdss6)));
subplot(1,8,7)
plot(stdss7)

stdss8 = stdss7(stdss7<(mean(stdss7,'all')+3*std(stdss7)));
subplot(1,8,8)
plot(stdss8)



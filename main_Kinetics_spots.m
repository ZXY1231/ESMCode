table_name = "D:\20211017\trackmate\5nmIgG_5nmIgA_4ms_200fps_1fold_avg10_diff_img_seletion_8_5_2_5\All Spots statistics.csv";
data1 = readtable(table_name);

data01(:,1) = data1(:,9); % Frame Number
data01(:,2) = data1(:,5); % x position
data01(:,3) = data1(:,6); % y position
data01(:,4) = data1(:,17); % Total intensity
data01(:,5) = data1(:,4); % Quality

data = table2array(data01);

data(:,1) = data(:,1)+1; % frame

% quality selection
% p = 0.3<data(:,5);
% data = data(p,:);
% q =  0.4>data(:,5);
% data = data(q,:);


%save by frames
all_spots = cell(1, data(end,1));
for i = 1:size(data,1)
    frame = data(i,1);
    all_spots{frame}(end+1,:) = data(i,:);
end

%% position diff<1 && time diff < 20 == the same particle
% active_trackers = TrackerInitializaon(all_spots{1},1);
active_trackers = {};
done_trackers = {};
new_trackers = {};
time_length = size(all_spots,2);

tic
for i = 1:time_length-1
    if isempty(all_spots{i})
        done_trackers = [done_trackers active_trackers]; 
        active_trackers = {};
        continue
    end
    
    if isempty(active_trackers)
        active_trackers = TrackerInitializaon(all_spots{i}(:,1:5),i);
        continue
    end
    
    spots = all_spots{i};
    
    spots = spots(:,1:5);% Frame,x,y,I,Q
    done_trackers_idx = [];
    used_spots = [];
    
    for j  = 1:size(active_trackers,2)
        tracker_xy = active_trackers{j}.positions(end,2:3);
        find_spot = active_trackers{j}.FindNext(tracker_xy,spots(:,2:3),2);
        if find_spot(1)
            active_trackers{j}.positions(end+1,:) = spots(find_spot(2),1:5); 
            used_spots(end+1) = find_spot(2);
        else
            done_trackers{end+1} = active_trackers{j};
            done_trackers_idx(end+1) = j;
        end
    end
    all_spots{i}(used_spots,:) = [];
    active_trackers(done_trackers_idx) = [];
    
    new_trackers = AllNewTracks(size(done_trackers,2)+size(active_trackers,2), all_spots{i}(:,1:5), i);
    active_trackers = [active_trackers new_trackers];     
end
%%
save('done_trackers_6_026')
toc

%% remove close spots in one frame

for i = 1:size(all_spots,2)
    one_frame_spots = all_spots{i};
    if size(one_frame_spots,1) <= 1
        continue
    end
    
    spots_xy = one_frame_spots(:,2:3);
    spots_dist = squareform(pdist(spots_xy));
    spots_dist = spots_dist+100*eye(size(spots_dist));
    [x,y] = find(spots_dist<30);
    one_frame_spots(x,:) = [];
    all_spots{i} = one_frame_spots;
end

%%
intensity = zeros(size(done_trackers,2),5);
for i = 1:size(done_trackers,2)
    [M, Idx] = max(done_trackers{i}.positions(:,4));
    intensity(i,:) = done_trackers{i}.positions(Idx,:);
end
% save('intensity_6_026')
% intensity2 = [];
% for i = 1:size(done_trackers,2)
%     if size(done_trackers{i}.positions,1)<15
%         continue
%     end
%     [M, Idx] = max(done_trackers{i}.positions(:,4));
%     intensity2(end+1,:) = done_trackers{i}.positions(Idx,:);
% end

% for i = 1:size(done_trackers,2)
%     one_particle = done_trackers{i};
%     
%     [M, Idx] = max(done_trackers{i}.positions(:,4));
%     intensity2(end+1,:) = done_trackers{i}.positions(Idx,:);
% end
%% remove bounday particles
boundary_particles = [];
for i = 1:size(intensity,1)
    x = intensity(i,2);
    y = intensity(i,3);
    if x<3||x>645||y<3||y>485
        boundary_particles(end+1) = i;
    end
end
intensity(boundary_particles,:) = [];
%%
source_path = 'D:\20210518\SLD135mA_5nmIgA_5nmantiIgA_05ms_4fold_smooth25avg_diff\';
imges = dir(append(source_path ,'*.tif'));
imge_num = length(imges);
shapes = size(imread(append(imges(1).folder, '/' ,imges(1).name)));
%%

roi = ones(6,6);
h = size(roi,1);
w = size(roi,2);
slip_win = 19;
all_images = zeros(slip_win, shapes(1), shapes(2),'single');
[xx, yy] = meshgrid(1:h, 1:w);
XY(:,:,1)=xx;
XY(:,:,2)=yy;
lb1 = [1,1,1,0.1,0.1];
ub1 = [inf,5,5,6,6];
lb2 = [0,1,0];
ub2 = [inf,slip_win,inf];
func1 = @(var,x) (var(1)*exp(-(x(:,:,1)-var(2)).^2/(2*var(4)^2)-(x(:,:,2)-var(3)).^2/(2*var(5)^2))); 
func2 = @(var,x) (-var(1)*abs(x-var(2))+var(3));  
options = optimset('Display','off');
% options = optimset('Display','off');
fit_intensity = intensity;
fit_intensity(:,4) = 0;
for i = 2:size(intensity,1)
    temp = num2cell(intensity(i,:));
    [mid_frame, x, y, inten] = temp{:};
    x = round(x);
    y = round(y);
    for j = 1:slip_win
        idx = mid_frame-floor(slip_win/2)+j-1;
        img = imread(append(imges(idx).folder, '/' ,imges(idx).name));
        all_images(j,:,:) = img;
    end
    patterns = all_images(:,y-h/2+1:y+h/2, x-w/2+1:x+w/2);
    inten_slip_win = zeros(1,slip_win);
    for j = 1:slip_win
        particle_pattern = squeeze(patterns(j,:,:));
        pattern_max = double(max(particle_pattern,[],'all'));
        ub1 = [pattern_max*2,5,5,6,6];
        result1 = lsqcurvefit(func1,...
            [pattern_max,h/2,w/2,h/6,w/6],XY,double(particle_pattern),...
            lb1,ub1,options);
        inten_slip_win(j) = result1(1);
    end
    result2 = lsqcurvefit(func2,[10,slip_win/2,10],1:slip_win,inten_slip_win,lb2,ub2,options);
    fit_intensity(i,4) = result2(3);
    if i<1&&i<10
        figure(1)
        scatter(1:slip_win,inten_slip_win)
        hold on
        plot(1:0.1:slip_win,func2(result2,1:0.1:slip_win))
        legend(int2str(i))
        hold off
        pause;   
    end

end
%%
save('fit_intensity_6_026')

%%
test = func2(result2,1:slip_win);
figure(12)
scatter(1:slip_win,inten_slip_win)
hold on
plot(1:0.1:slip_win,func2(result2,1:0.1:slip_win))
toc

%%
data = intensity;
time_window = 30000;
figure(13)
for i = 1:15
    subplot(3,5,i)
    p = ((i-1)*time_window<data(:,1))&(data(:,1)<i*time_window);
    data_window = data(p,:);
    histogram(data_window(:,4),-5:2:80)
%     histogram(data_window(:,4))
    legend(int2str(size(data_window(:,4),1)))
%     axis([-5 40 0 500])
end

%%
data = fit_intensity;
time_window = 30000;
figure(14)
for i = 1:15
    subplot(3,5,i)
    p = ((i-1)*time_window<data(:,1))&(data(:,1)<i*time_window);
    data_window = data(p,:);
    histogram(data_window(:,4),-5:1:20)
%     histogram(data_window(:,4))
    legend(int2str(size(data_window(:,4),1)))
%     axis([-5 40 0 500])
end

%%
figure(15)
histogram(fit_intensity(:,4),0:0.1:15)
legend(int2str(size(data(:,4),1)))
%%
% a = data(:,1)<1125;
% test1 = data(a,:);
% time_window = 4*375;
% figure(4)
% for i = 1:20
%     subplot(4,5,i)
%     p = ((i-1)*time_window<data(:,1))&(data(:,1)<i*time_window);
%     data_window = data(p,:);
%     histogram(data_window(:,4),0:2:80)
% end
%%
figure(18)
scatter(data(:,4),data(:,1))
% axis([0 80 0 20])
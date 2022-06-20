 tic
trackmate_tracks = [ "C:\Users\xzhou145\Desktop\20220309\trackmate\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_000001_070000_8_1_5_2_2_1\Spots in tracks statistics.csv",...
                    ];
                
%                 C:\Users\xzhou145\Desktop\20220309\trackmate\chip1_2_5nMIgGA_6ms_160fps_1fold_smoothavg24_diff_004000_084000_8_1_5_2_2_1\Spots in tracks statistics.csv

particle_info = [];
frame_min = 9600;
frame_sec = 160;
duration = 70001;
noise_frame_spots_thre = 60;
track_range = [24 500];
% track_range = [8 500];
close_dis_thre = 14;
t_window = 60; % s
t_step = 5;
monomer_con = 2.5; %nM
view_range = [5 128 25 75];
% close_spots_thre = 6;

%%
for folder = 1:1
    folder
    clear tracks_data1 tracks_data01
    table_name = trackmate_tracks(folder);
    tracks_data1 = readtable(table_name);

    tracks_data01(:,1) = tracks_data1(:,9); % Frame Number
    tracks_data01(:,2) = tracks_data1(:,3); % track ID
    tracks_data01(:,3) = tracks_data1(:,5); % x position
    tracks_data01(:,4) = tracks_data1(:,6); % y position
    tracks_data01(:,5) = tracks_data1(:,17); % Total intensity
    tracks_data01(:,6) = tracks_data1(:,4); % Quality
    tracks_data01(:,7) = tracks_data1(:,16); % Max_Intensity

    tracks_data = table2array(tracks_data01);

    tracks_data(:,1) = tracks_data(:,1)+1; % frame
    tracks_data(:,2) = tracks_data(:,2)+1; % track ID

    %save by tracks
    all_tracks = cell(1,tracks_data(end,2));
    for i = 1:size(tracks_data,1)
        id_idx = tracks_data(i,2);
        all_tracks{id_idx}(end+1,:) = tracks_data(i,:); 
    end
    
    %save by frames
    all_frames = cell(1, max(tracks_data(:,1)));
    for i = 1:size(tracks_data,1)
        frame = tracks_data(i,1);
        all_frames{frame}(end+1,:) = tracks_data(i,:); % save IDs by frame
    end
    
    %% remove abnormal frame associdated tracks, too many noisy spots
    abnormal_tracks = ones(size(all_tracks));
    for i = 1:size(all_frames,2)
        frame_info = all_frames{i};
        if size(frame_info,1)>noise_frame_spots_thre
            IDs = frame_info(:,2);
            abnormal_tracks(IDs) = 0; 
        end
    end
    abnormal_frames = [];
    for i = 1:size(abnormal_tracks,2)
        if abnormal_tracks(1,i) == 1
            continue
        end
        one_track = all_tracks{i};
        abnormal_frames = cat(2, abnormal_frames, one_track(:,1)');
    end
    abnormal_frames = unique(abnormal_frames);
    for i = 1:size(abnormal_frames,2)
        if size(abnormal_frames,1) == 0
            continue
        end
        frame_info = all_frames{abnormal_frames(1,i)};
        IDs = frame_info(:,2);
        abnormal_tracks(IDs) = 0; 
    end
%     all_tracks = all_tracks(logical(abnormal_tracks));
    
    %% remove close spots in one frame
% method1 
%     for i = 1:size(all_frames,2)
%         one_frame_spots = all_frames{i};
%         if size(one_frame_spots,1) <= 1
%             continue
%         end
%         spots_xy = one_frame_spots(:,2:3);
%         spots_dist = squareform(pdist(spots_xy));
%         spots_dist = spots_dist+100*eye(size(spots_dist));
%         [x,y] = find(spots_dist<close_dis_thre);
%         close_spot_x = unique(x);
%         if length(close_spot_x) > close_spots_thre
%             one_frame_spots(x,:) = [];
%             all_spots{i} = one_frame_spots;
%         else
%             max_supress = 1;
%             for k = 2:length(close_spot_x)/2      % x here is the x axis of full metrix
%                 if one_frame_spots(k,4)>one_frame_spots(max_supress,4)
%                     max_supress = k;
%                 end
%             end
%             one_frame_spots = one_frame_spots(max_supress,:);
%             all_spots{i} = one_frame_spots;
%         end
%     end
% method2
%     for i = 1:size(all_frames,2)
%     for i = 1:3896
%         one_frame_spots = all_frames{i};
%         if size(one_frame_spots,1) <= 1
%             continue
%         end
%         spots_xy = one_frame_spots(:,2:3);
%         spots_dist = squareform(pdist(spots_xy));
%         spots_dist = spots_dist+100*eye(size(spots_dist));
%         [x,y] = find(spots_dist<close_dis_thre);
%         close_spot_x = unique(x);
%         if length(close_spot_x) <=1 %% no closed spots
%             continue
%         elseif length(close_spot_x) > close_spots_thre %% extensive closed spots
%             IDs = one_frame_spots(close_spot_x,2);
%             abnormal_tracks(IDs) = 0; 
%         else
%             max_supress = close_spot_x(1);
%             for k = 2:length(close_spot_x)      
%                 if one_frame_spots(k,4)>one_frame_spots(max_supress,4)
%                     max_supress = k;
%                 end
%             end
%             one_frame_spots = one_frame_spots(max_supress,:);
%             all_spots{i} = one_frame_spots;
%         end
%     end

% need to be improved
    close_frames_check = zeros(size(all_frames));
    one_group_close_IDs = [];
%     size(all_tracks,2)
    for i = 1:size(all_tracks,2)
        one_track = all_tracks{i};
        one_track_frames = one_track(:,1);
        
        for j = 1:length(one_track_frames)
            if close_frames_check(one_track_frames(j)) == 1
                continue
            end
            one_frame_spots = all_frames{one_track_frames(j)};
            if size(one_frame_spots,1) <= 1
                continue
            end
            spots_xy = one_frame_spots(:,3:4);
            spots_dist = squareform(pdist(spots_xy));
            spots_dist = spots_dist+100*eye(size(spots_dist));
            [x,~] = find(spots_dist<close_dis_thre);
            close_spot_x = unique(x);
            if length(close_spot_x) <=1 %% no closed spots
                continue
            else 
                IDs = one_frame_spots(close_spot_x,2);
                one_group_close_IDs = cat(1,one_group_close_IDs,IDs);
            end
            close_frames_check(one_track_frames(j)) = 1;
        end
        
        one_group_close_IDs = unique(one_group_close_IDs);
        if length(one_group_close_IDs)<=1
            continue
        end
        
        max_supress_ID = one_group_close_IDs(1);
        max_supress = 1;
        for k = 2:length(one_group_close_IDs)      
            if max(all_tracks{one_group_close_IDs(k)}(:,5))>...
                    max(all_tracks{max_supress_ID}(:,5))
                max_supress_ID = one_group_close_IDs(k);
                max_supress = k;
            end
        end
        one_group_close_IDs(max_supress) = [];
        abnormal_tracks(one_group_close_IDs) = 0;
        one_group_close_IDs = [];
    end
    all_tracks = all_tracks(logical(abnormal_tracks));
%%

    %% track length band pass
    intensity = zeros(size(all_tracks,2),7);
    for i = 1:size(all_tracks,2)
        if size(all_tracks{i},1)<track_range(1)        
            continue

        end
        if size(all_tracks{i},1)>track_range(2)
            continue

        end
        [M, Idx] = max(all_tracks{i}(:,5));
        intensity(i,:) = all_tracks{i}(Idx,:);
    end
    intensity(intensity(:,1)==0,:) = [];
    %% remove bounday particles
    boundary_particles = [];
    for i = 1:size(intensity,1)
        x = intensity(i,3);
        y = intensity(i,4);
        if x<view_range(1)||x>view_range(2)||y<view_range(3)||y>view_range(4)
            boundary_particles(end+1) = i;
        end
    end
    intensity(boundary_particles,:) = [];
%     
    intensity(:,1) = intensity(:,1)+(folder-1)*frame_min; % correct frames
    particle_info = cat(1,particle_info, intensity);
    toc
end

%% plot intensity counts f(t)
range = [2, 2.2, 2.45, 3.2, 5, 5, 100];
range = [2, 2.8, 2.9, 5.1, 8, 8, 100];
range = [2, 2.2, 2.45, 3.25, 8, 8, 100];
range = [1, 2, 3.8, 3.8, 5.8, 5.8, 100];
range = [1, 1.5, 2.4, 3.0, 4.6, 4.6, 100];
range = [1, 1.5, 2.4, 3.0, 4.6, 4.6, 5.8,5.8,100];
range = [1, 1.5, 3.2, 3.2, 4.6, 4.6, 6.8,6.8,100];
% range = [1, 1.5, 3.2, 3.2, 4.6, 4.6, 10.2,10.2,100];

counts = zeros(4,duration);

for i = 1:size(particle_info,1)
    particle = particle_info(i,:);
    qua = particle(1,6);
    if qua>range(9)
        continue
    elseif qua>=range(8)
        counts(4,particle(1)) = counts(4,particle(1,1))+1;
    elseif qua>=range(6) && qua<range(7)
        counts(3,particle(1)) = counts(3,particle(1,1))+1;
    elseif qua>=range(4) && qua<range(5)
        counts(2,particle(1)) = counts(2,particle(1,1))+1;
    elseif qua>=range(2) && qua<range(3)
        counts(1,particle(1)) = counts(1,particle(1,1))+1;
    end
end


%% total counts total intensity
q_cut = range;
p_slcd = (particle_info(:,6)>=q_cut(2))&(particle_info(:,6)<q_cut(3));
particle_info_s = particle_info(p_slcd,:);
p_slcd = (particle_info(:,6)>=q_cut(4))&(particle_info(:,6)<=q_cut(5));
particle_info_s2 = particle_info(p_slcd,:);
p_slcd = (particle_info(:,6)>=q_cut(6))&(particle_info(:,6)<=q_cut(7));
particle_info_s3 = particle_info(p_slcd,:);
p_slcd = (particle_info(:,6)>=q_cut(8))&(particle_info(:,6)<=q_cut(9));
particle_info_s4 = particle_info(p_slcd,:);

figure(11)
subplot(1,5,1)
histogram(particle_info(:,5),-400:35:3600)
legend(int2str(size(particle_info(:,5),1)))
title('total')
subplot(1,5,2)
histogram(particle_info_s(:,5),-400:40:1400)
legend(int2str(size(particle_info_s(:,5),1)))
title('IgG')
subplot(1,5,3)
histogram(particle_info_s2(:,5),-400:116:2200)
legend(int2str(size(particle_info_s2(:,5),1)))
title('IgA')
subplot(1,5,4)
histogram(particle_info_s3(:,5),-400:40:3600)
legend(int2str(size(particle_info_s3(:,5),1)))
title('Complex')
subplot(1,5,5)
histogram(particle_info_s4(:,5),-400:80:3600)
legend(int2str(size(particle_info_s4(:,5),1)))
title('IgM')
%% total counts max intensity
q_cut = range;
p_slcd = (particle_info(:,6)>=q_cut(2))&(particle_info(:,6)<q_cut(3));
particle_info_s = particle_info(p_slcd,:);
p_slcd = (particle_info(:,6)>=q_cut(4))&(particle_info(:,6)<=q_cut(5));
particle_info_s2 = particle_info(p_slcd,:);
p_slcd = (particle_info(:,6)>=q_cut(6))&(particle_info(:,6)<=q_cut(7));
particle_info_s3 = particle_info(p_slcd,:);
p_slcd = (particle_info(:,6)>=q_cut(8))&(particle_info(:,6)<=q_cut(9));
particle_info_s4 = particle_info(p_slcd,:);
figure(1)
subplot(1,5,1)
histogram(particle_info(:,7),-40:1.5:240)
legend(int2str(size(particle_info(:,7),1)))
title('total')
subplot(1,5,2)
histogram(particle_info_s(:,7),-40:4:80)
legend(int2str(size(particle_info_s(:,7),1)))
title('IgG')
subplot(1,5,3)
histogram(particle_info_s2(:,7),-40:4:120)
legend(int2str(size(particle_info_s2(:,7),1)))
title('IgA')
subplot(1,5,4)
histogram(particle_info_s3(:,7),-40:4:200)
legend(int2str(size(particle_info_s3(:,5),1)))
title('Complex')
subplot(1,5,5)
histogram(particle_info_s4(:,7),-40:4:220)
legend(int2str(size(particle_info_s4(:,5),1)))
title('IgM')

%% counts accmulate_window
counts_cum = cumsum(counts,2);

counts_sec = 60*[1:duration]/frame_min;

figure(13)
hold on
plot(counts_sec,counts_cum(1,:))
plot(counts_sec,counts_cum(2,:))
plot(counts_sec,counts_cum(3,:))
plot(counts_sec,counts_cum(4,:))
hold off
title('counts accumulation')
xlabel('time /s')
ylabel('counts')
% axis([1 600 0 inf])
legend('IgG','IgA','Complex','IgM')

%% time window

counts_window = zeros(4,floor((duration-frame_sec*t_window)/(frame_sec*t_step)));
% counts_window = zeros(4,floor((duration/frame_min)*(60/t_window)));

for i = 1:size(counts_window,2)
    t = frame_sec*(i-1)*t_step+1:frame_sec*((i-1)*t_step+t_window);
    counts_window(:,i) = sum(counts(:,t),2);
end


figure(14)
hold on
plot(t_window:t_step:t_step*(size(counts_window,2)-1)+t_window,counts_window(1,:))
plot(t_window:t_step:t_step*(size(counts_window,2)-1)+t_window,counts_window(2,:))
plot(t_window:t_step:t_step*(size(counts_window,2)-1)+t_window,counts_window(3,:))
plot(t_window:t_step:t_step*(size(counts_window,2)-1)+t_window,counts_window(4,:))
hold off
title('counts along the time')
xlabel('time /min')
ylabel('counts')
% axis([0 duration-1 0 inf])
% axis([1.5 duration-1 0 inf])
legend('IgG', 'IgA', 'Complex','IgM','Location','east')

%%
% Kd = zeros(1,20);
% f_conv_all = zeros(1,20);
% for t = 1:20
%     f_conv = 5/(counts_window(2,t)+counts_window(3,t));
%     f_conv_all(t) = f_conv;
%     IgA_b = f_conv*counts_window(3,t);
%     IgA_unb = f_conv*counts_window(2,t);
%     IgG = 5-IgA_b;
%     Kd(t) = IgA_unb*IgG/IgA_b;
% end
% figure(103)
% plot(0.5:0.5:10,Kd)
% axis([0 10 0 inf])
%% k fitting 2
antiIgA_counts = counts_window(1,:);
IgA_counts = counts_window(2,:);
Complex_counts = counts_window(3,:);
ratio = IgA_counts./(IgA_counts+Complex_counts);

ratio = (IgA_counts./(IgA_counts+Complex_counts)+antiIgA_counts./(antiIgA_counts+Complex_counts))/2;
% t = 90:30:960;
t = t_window:t_step:t_step*(size(counts_window,2)-1)+t_window;
t = t+50;
% t = t();
ratio = ratio(1:length(t));
func = @(var,x) (((1+var(2))*(1-var(3))./(1+var(2)*exp(x/var(1))))+var(3)); %exp ratio change
options = optimset('Display','off');
result1 = lsqcurvefit(func,...
    [50,1,0.5],t,ratio,...
    [],[],options);

% use nM unit here
func2 = @(var,x)((1+var(2))*(var(3)-1)*(var(2)/var(1))*exp(x/var(1))./((1+var(2)*exp(x/(var(1)))).^2));%derivative of func
y = func2(result1,t)*monomer_con;
func3 = @(var,x) (var(1)*(monomer_con-x)-var(2)*(x).^2);
IgA_con = ratio*monomer_con;
result2 = lsqcurvefit(func3,...
    [1e-2,1e-2],IgA_con,y,...
    [0 0],[],options)
kd = result2(1)/result2(2)

%%
t2 = t(1):1:t(end);
y2 = func(result1,t2);
figure(220)
scatter(t-t_window/2,ratio)
hold on
plot(t2-t_window/2,y2)
hold off
xlabel('time /s')
% ylabel('[antiIgA-t]/[antiIgA-total]')
ylabel('[antiIgA-t]/[antiIgA-total]')
ylabel('[IgA-t]/[IgA-total]')
ylabel('Average of antiIgA and IgA ratio')
% axis([50 300 0.65 1])

tic
% trackmate_tracks_bind = [ "C:\Users\xzhou145\Desktop\20220309\trackmate\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_080000_144000_8_1_5_2_2_1\Spots in tracks statistics.csv",...
%                     ];
% 
% trackmate_tracks_unbind = [ "C:\Users\xzhou145\Desktop\20220309\trackmate\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_080000_144000_pixel_inverse_8_1_5_2_2_1\Spots in tracks statistics.csv",...
%                     ];
                
                
trackmate_tracks_bind = [ "C:\Users\xzhou145\Desktop\20220309\trackmate\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_000001_070000_8_1_5_2_2_1\Spots in tracks statistics.csv",...
                    ];

trackmate_tracks_unbind = [ "C:\Users\xzhou145\Desktop\20220309\trackmate\chip2_2_5nMIgM_6ms_160fps_1fold_smoothavg24_diff_000001_070000_pixel_inverse_8_1_5_2_2_1\Spots in tracks statistics.csv",...
                    ];

particle_info_bind = [];
particle_info_unbind = [];
frame_min = 9600;
frame_sec = 160;
duration = 70001;
noise_frame_spots_thre = 60;
track_range = [24 50];
close_dis_thre = 14;
t_window = 60; % s
t_step = 10;
monomer_con = 2.5; %nM
bind_time_max = 50;
bind_unbind_dis = 16;
view_range = [10 120 25 75];
view_range = [5 128 5 75];

%% binding particles
for folder = 1:1
    folder
    clear tracks_data1 tracks_data01
    table_name_bind = trackmate_tracks_bind(folder);
    tracks_data1_bind = readtable(table_name_bind);

    tracks_data01_bind(:,1) = tracks_data1_bind(:,9); % Frame Number
    tracks_data01_bind(:,2) = tracks_data1_bind(:,3); % track ID
    tracks_data01_bind(:,3) = tracks_data1_bind(:,5); % x position
    tracks_data01_bind(:,4) = tracks_data1_bind(:,6); % y position
    tracks_data01_bind(:,5) = tracks_data1_bind(:,17); % Total intensity
    tracks_data01_bind(:,6) = tracks_data1_bind(:,4); % Quality
    tracks_data01_bind(:,7) = tracks_data1_bind(:,16); % Max_Intensity

    tracks_data_bind = table2array(tracks_data01_bind);

    tracks_data_bind(:,1) = tracks_data_bind(:,1)+1; % frame
    tracks_data_bind(:,2) = tracks_data_bind(:,2)+1; % track ID

    %save by tracks
    all_tracks_bind = cell(1,tracks_data_bind(end,2));
    for i = 1:size(tracks_data_bind,1)
        id_idx_bind = tracks_data_bind(i,2);
        all_tracks_bind{id_idx_bind}(end+1,:) = tracks_data_bind(i,:); 
    end
    
    %save by frames
    all_frames_bind = cell(1, max(tracks_data_bind(:,1)));
    for i = 1:size(tracks_data_bind,1)
        frame = tracks_data_bind(i,1);
        all_frames_bind{frame}(end+1,:) = tracks_data_bind(i,:); % save IDs by frame
    end
    
    %% remove abnormal frame associdated tracks
    abnormal_tracks_bind = ones(size(all_tracks_bind));
    for i = 1:size(all_frames_bind,2)
        frame_info_bind = all_frames_bind{i};
        if size(frame_info_bind,1)>noise_frame_spots_thre
            IDs = frame_info_bind(:,2);
            abnormal_tracks_bind(IDs) = 0; 
        end
    end
    abnormal_frames_bind = [];
    for i = 1:size(abnormal_tracks_bind,2)
        if abnormal_tracks_bind(1,i) == 1
            continue
        end
        one_track = all_tracks_bind{i};
        abnormal_frames_bind = cat(2, abnormal_frames_bind, one_track(:,1)');
    end
    abnormal_frames_bind = unique(abnormal_frames_bind);
    for i = 1:size(abnormal_frames_bind,2)
        if size(abnormal_frames_bind,1) == 0
            continue
        end
        frame_info_bind = all_frames_bind{abnormal_frames_bind(1,i)};
        IDs = frame_info_bind(:,2);
        abnormal_tracks_bind(IDs) = 0; 
    end
%     all_tracks_bind = all_tracks_bind(logical(abnormal_tracks_bind));
    
    %% remove close spots in one frame

% need to be improved
    close_frames_check = zeros(size(all_frames_bind));
    one_group_close_IDs = [];
%     size(all_tracks,2)
    for i = 1:size(all_tracks_bind,2)
        one_track = all_tracks_bind{i};
        one_track_frames = one_track(:,1);
        
        for j = 1:length(one_track_frames)
            if close_frames_check(one_track_frames(j)) == 1
                continue
            end
            one_frame_spots = all_frames_bind{one_track_frames(j)};
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
            if max(all_tracks_bind{one_group_close_IDs(k)}(:,5))>...
                    max(all_tracks_bind{max_supress_ID}(:,5))
                max_supress_ID = one_group_close_IDs(k);
                max_supress = k;
            end
        end
        one_group_close_IDs(max_supress) = [];
        abnormal_tracks_bind(one_group_close_IDs) = 0;
        one_group_close_IDs = [];
    end
    all_tracks_bind = all_tracks_bind(logical(abnormal_tracks_bind));
%%

    %% track length bang pass
    intensity_bind = zeros(size(all_tracks_bind,2),7);
    for i = 1:size(all_tracks_bind,2)
        if size(all_tracks_bind{i},1)<track_range(1)
            continue
        end
        if size(all_tracks_bind{i},1)>track_range(2)
            continue
        end
        [M, Idx] = max(all_tracks_bind{i}(:,5));
        intensity_bind(i,:) = all_tracks_bind{i}(Idx,:);
    end
    intensity_bind(intensity_bind(:,1)==0,:) = [];
    %% remove bounday particles
    boundary_particles_bind = [];
    for i = 1:size(intensity_bind,1)
        x = intensity_bind(i,3);
        y = intensity_bind(i,4);
        if x<view_range(1)||x>view_range(2)||y<view_range(3)||y>view_range(4)
            boundary_particles_bind(end+1) = i;
        end
    end
    intensity_bind(boundary_particles_bind,:) = [];
%     
    intensity_bind(:,1) = intensity_bind(:,1)+(folder-1)*frame_min; % correct frames
    particle_info_bind = cat(1,particle_info_bind, intensity_bind);
    toc
end

%% unbinding particles
for folder = 1:1
    folder
    clear tracks_data1 tracks_data01
    table_name_unbind = trackmate_tracks_unbind(folder);
    tracks_data1_unbind = readtable(table_name_unbind);

    tracks_data01_unbind(:,1) = tracks_data1_unbind(:,9); % Frame Number
    tracks_data01_unbind(:,2) = tracks_data1_unbind(:,3); % track ID
    tracks_data01_unbind(:,3) = tracks_data1_unbind(:,5); % x position
    tracks_data01_unbind(:,4) = tracks_data1_unbind(:,6); % y position
    tracks_data01_unbind(:,5) = tracks_data1_unbind(:,17); % Total intensity
    tracks_data01_unbind(:,6) = tracks_data1_unbind(:,4); % Quality
    tracks_data01_unbind(:,7) = tracks_data1_unbind(:,16); % Max_Intensity

    tracks_data_unbind = table2array(tracks_data01_unbind);

    tracks_data_unbind(:,1) = tracks_data_unbind(:,1)+1; % frame
    tracks_data_unbind(:,2) = tracks_data_unbind(:,2)+1; % track ID

    %save by tracks
    all_tracks_unbind = cell(1,tracks_data_unbind(end,2));
    for i = 1:size(tracks_data_unbind,1)
        id_idx_unbind = tracks_data_unbind(i,2);
        all_tracks_unbind{id_idx_unbind}(end+1,:) = tracks_data_unbind(i,:); 
    end
    
    %save by frames
    all_frames_unbind = cell(1, max(tracks_data_unbind(:,1)));
    for i = 1:size(tracks_data_unbind,1)
        frame = tracks_data_unbind(i,1);
        all_frames_unbind{frame}(end+1,:) = tracks_data_unbind(i,:); % save IDs by frame
    end
    
    %% remove abnormal frame associdated tracks
    abnormal_tracks_unbind = ones(size(all_tracks_unbind));
    for i = 1:size(all_frames_unbind,2)
        frame_info_unbind = all_frames_unbind{i};
        if size(frame_info_unbind,1)>noise_frame_spots_thre
            IDs = frame_info_unbind(:,2);
            abnormal_tracks_unbind(IDs) = 0; 
        end
    end
    abnormal_frames_unbind = [];
    for i = 1:size(abnormal_tracks_unbind,2)
        if abnormal_tracks_unbind(1,i) == 1
            continue
        end
        one_track = all_tracks_unbind{i};
        abnormal_frames_unbind = cat(2, abnormal_frames_unbind, one_track(:,1)');
    end
    abnormal_frames_unbind = unique(abnormal_frames_unbind);
    for i = 1:size(abnormal_frames_unbind,2)
        if size(abnormal_frames_unbind,1) == 0
            continue
        end
        frame_info_unbind = all_frames_unbind{abnormal_frames_unbind(1,i)};
        IDs = frame_info_unbind(:,2);
        abnormal_tracks_unbind(IDs) = 0; 
    end
%     all_tracks_unbind = all_tracks_unbind(logical(abnormal_tracks_unbind));
    
    %% remove close spots in one frame

% need to be improved
    close_frames_check = zeros(size(all_frames_unbind));
    one_group_close_IDs = [];
%     size(all_tracks,2)
    for i = 1:size(all_tracks_unbind,2)
        one_track = all_tracks_unbind{i};
        one_track_frames = one_track(:,1);
        
        for j = 1:length(one_track_frames)
            if close_frames_check(one_track_frames(j)) == 1
                continue
            end
            one_frame_spots = all_frames_unbind{one_track_frames(j)};
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
            if max(all_tracks_unbind{one_group_close_IDs(k)}(:,5))>...
                    max(all_tracks_unbind{max_supress_ID}(:,5))
                max_supress_ID = one_group_close_IDs(k);
                max_supress = k;
            end
        end
        one_group_close_IDs(max_supress) = [];
        abnormal_tracks_unbind(one_group_close_IDs) = 0;
        one_group_close_IDs = [];
    end
    all_tracks_unbind = all_tracks_unbind(logical(abnormal_tracks_unbind));
%%

    %% track length bang pass
    intensity_unbind = zeros(size(all_tracks_unbind,2),7);
    for i = 1:size(all_tracks_unbind,2)
        if size(all_tracks_unbind{i},1)<track_range(1)
            continue
        end
        if size(all_tracks_unbind{i},1)>track_range(2)
            continue
        end
        [M, Idx] = max(all_tracks_unbind{i}(:,5));
        intensity_unbind(i,:) = all_tracks_unbind{i}(Idx,:);
    end
    intensity_unbind(intensity_unbind(:,1)==0,:) = [];
    %% remove bounday particles
    boundary_particles_unbind = [];
    for i = 1:size(intensity_unbind,1)
        x = intensity_unbind(i,3);
        y = intensity_unbind(i,4);
        if x<view_range(1)||x>view_range(2)||y<view_range(3)||y>view_range(4)
            boundary_particles_unbind(end+1) = i;
        end
    end
    intensity_unbind(boundary_particles_unbind,:) = [];
%     
    intensity_unbind(:,1) = intensity_unbind(:,1)+(folder-1)*frame_min; % correct frames
    particle_info_unbind = cat(1,particle_info_unbind, intensity_unbind);
    toc
end

%% find transient event, code efficiency need to be improved
frames_unbind = particle_info_unbind(:,1);
particle_transient_bind = zeros(size(particle_info_bind,1),1);
particle_transient_unbind = zeros(size(particle_info_unbind,1),1);
for i = 1:size(particle_info_bind,1)
    particle_bind = particle_info_bind(i,:);
    frame = particle_bind(1);
    frames_cand = find((frames_unbind-frame)<bind_time_max&(frames_unbind-frame)>0);
    particles_unbind_cand = particle_info_unbind(frames_cand,:);
    for j = 1:size(particles_unbind_cand,1)
        p_u_cand = particles_unbind_cand(j,:);
        if ((particle_bind(3)-p_u_cand(3))^2+(particle_bind(4)-p_u_cand(4))^2)<bind_unbind_dis
            particle_transient_unbind(frames_cand(j)) = 1;
            particle_transient_bind(i) = 1;
        end
    end
    
end

transient_particles_unbind = particle_info_unbind(logical(particle_transient_unbind),:);
transient_particles_bind = particle_info_bind(logical(particle_transient_bind),:);

nontransient_particles_unbind = particle_info_unbind(logical(-particle_transient_unbind+1),:);
nontransient_particles_bind = particle_info_bind(logical(-particle_transient_bind+1),:);

%% plot intensity counts f(t)
particle_info = nontransient_particles_unbind;
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

figure(3)
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

% %% delete transient event in the IgG IgA counts, run main_Kinetics_multifolder_trackmatetracks2 firstly
% [val,pos]=intersect(transient_particles_bind(:,1),particle_info_s(:,1));
% [val2,pos2]=intersect(transient_particles_bind(:,1),particle_info_s2(:,1));
% particle_info_s(pos,:) = [];
% particle_info_s2(pos2,:) = [];
% 
% counts(1,val) = counts(1,val)-1;
% counts(2,val) = counts(2,val)-1;

% %% add transient event to complex counts
% for i = 1:size(transient_particles_bind,1)
%     particle = transient_particles_bind(i,:);
%     qua = particle(1,6);
%     counts(3,particle(1)) = counts(3,particle(1,1))+1;
% end

%%

%% time window

counts_window = zeros(4,floor((duration-frame_sec*t_window)/(frame_sec*t_step)));
% counts_window = zeros(4,floor((duration/frame_min)*(60/t_window)));

for i = 1:size(counts_window,2)
    t = frame_sec*(i-1)*t_step+1:frame_sec*((i-1)*t_step+t_window);
    counts_window(:,i) = sum(counts(:,t),2);
end


figure(4)
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
%% k fitting 2
IgA_counts = counts_window(2,:);
Complex_counts = counts_window(3,:);
ratio = IgA_counts./(IgA_counts+Complex_counts);

t = 90:30:960;
ratio = ratio(1:length(t));
func = @(var,x) (((1+var(2))*(1-var(3))./(1+var(2)*exp(x/var(1))))+var(3)); %exp ratio change
options = optimset('Display','off');
result1 = lsqcurvefit(func,...
    [50,1,0.5],t,ratio,...
    [],[],options);

% use nM unit here
func2 = @(var,x)((1+var(2))*(var(3)-1)*(var(2)/var(1))*exp(x/var(1))./((1+var(2)*exp(x/(var(1)))).^2));%derivative of func
y = func2(result1,t)*10;
func3 = @(var,x) (var(1)*(10-x)-var(2)*(x).^2);
IgA_con = ratio*10;
result2 = lsqcurvefit(func3,...
    [1e-2,1e-2],IgA_con,y,...
    [],[],options)
kd = result2(1)/result2(2)

%%
t2 = t(1)-60:1:t(end);
y2 = func(result1,t2);
figure(220)
scatter(t,ratio)
hold on
plot(t2,y2)
hold off
xlabel('time /s')
% ylabel('[antiIgA-t]/[antiIgA-total]')
ylabel('[IgA-t]/[IgA-total]')











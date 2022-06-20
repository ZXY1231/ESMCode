%Initialize trackers based on the first_frame_bright_spot
function all_tracks = TrackerInitializaon(spots,t)
    all_tracks = cell(1,size(spots,1));
    for i  = 1:size(spots,1)
        all_tracks{i} = Tracker(i,spots(i,:),0,t);
    end
end
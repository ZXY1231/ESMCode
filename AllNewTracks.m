function new_tracks = AllNewTracks(track_id_end, next_spots,t)

    new_tracks = cell(1,size(next_spots,1));
    for i = 1:size(next_spots,1)
        new_tracks{i} = Tracker(track_id_end+i, next_spots(i,:), 0, t);
    end
end
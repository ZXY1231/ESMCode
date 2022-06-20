function find_spot = TrackFindNext(one_tracker, next_spots)
    
    find_spot = one_tracker.FindNext(next_bright_particles);
    
    if find_spot(1) && all_images_bright_particles{t}(find_spot(2),3) == 0% multi-particles assignment temporary solution
%     if find_spot(1)
        all_images_bright_particles{t}(find_spot(2),3) = one_tracker.track_id;
    else 
        if find_spot(1)
            one_tracker.position_xy(end,:) = [];% delete multi-assigned pisitions
            
            
        end
        find_spot = one_tracker.FindNextDim(next_dim_particles);
        if find_spot(1)
            all_images_dim_particles{t}(find_spot(2),3) = one_tracker.track_id;
        end
    end
end
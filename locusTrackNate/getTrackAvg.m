function track_avg = getTrackAvg(tracks)
    
    ss = size(tracks);
        
    for itrack = 1:ss(1)
            
            track_sum = sum(tracks(itrack,~isnan(tracks(itrack,:))));
            num = sum(~isnan(tracks(itrack,:)));
            
            track_avg(itrack) = track_sum/num;
            
    end
    
end
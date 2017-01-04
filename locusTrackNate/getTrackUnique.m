function [track] = getTrackUnique(track1, track2)

%% THIS FUNTION FINDS UNIQUE VALUES BETWEEN TWO TRAJECOTIES, USED TO FIND
% THE SPLITTING POINT 
%
% NJK

ssL = size(track1);

epsil = 0.0001;


for iuniq = 1:ssL(2)
    
    if (track1(iuniq) - epsil) <= track2(iuniq) && track2(iuniq) <= (track1(iuniq) + epsil)
        
        track2(iuniq) = NaN;
        
    end
    
end

track = track2;

end
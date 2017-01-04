function [rel, flag] = getRelativeCoords( data ) 

global CONST

if isempty( CONST )
    if exist('loadConstantsMine','file');
        loadConstantsMine
    else
        loadConstants
    end
end


TimeStep     = CONST.getLocusTracks.TimeStep/60;
PixelSize    = CONST.getLocusTracks.PixelSize;



% FUNCTION TO GET CELL-RELATIVE COORDINATES

[L,L0] = getCellLength( data );

filename = ['Cell ', num2str(data.ID)];


if ~isfield(data,'track1') || ~isfield(data,'track2')
        
        disp([filename, ' : Getting New Locus Tracks']);
        tracks = getLocusTracksDev_lite(data);
        
    else
        
        disp([filename, ' : Using Existing Locus Tracks']);
        tracks = data;
        
    end
    
    % INITIALIZE A FEW VECTORS FOR GOOD MEASURE
    
    sort.green = [];
    sort.red = [];
    
    % SORT TRACKS IN XTRACK INTO UPPER AND LOWER TRACKS, DETERMINED BY THE
    % SIGN OF THE LAST 5 POSITION COORDINATES (NEGATIVE ==> LOWER)
    
    sort.green  = getTrackSort(tracks.track1.xtrack.*PixelSize);
    sort.red = getTrackSort(tracks.track2.xtrack.*PixelSize);
    
    % MAKE SURE THAT THERE ARE TWO OF EACH TRACK
    
    if isempty(sort.green.low) || isempty(sort.green.high)
        
        disp('Error: missing green track')
        
        g_init_rel = NaN;
        r_init_rel = NaN;
        g_high_rel = NaN;
        g_low_rel  = NaN;
        r_high_rel = NaN;
        r_low_rel  = NaN;
        
        flag = 1;
        
    else
        
        if isempty(sort.red.low) || isempty(sort.red.high)
            
            disp('Error: missing red track')
            
            g_init_rel = NaN;
            r_init_rel = NaN;
            g_high_rel = NaN;
            g_low_rel  = NaN;
            r_high_rel = NaN;
            r_low_rel  = NaN;
            
            flag = 1;
            
        else
            
            flag = 0;
            % MOST TRACKS ARE THE SAME, SO WILL GET ONE AVERAGE TRACK PER UPPER AND
            % LOWER ROUTE
            
            avg.green.low = getTrackAvg(sort.green.low);
            avg.green.high_ = getTrackAvg(sort.green.high);
            avg.red.low = getTrackAvg(sort.red.low);
            avg.red.high_ = getTrackAvg(sort.red.high);
            
            % REMOVE REPETITIVE VALUES IN TRAJECTORIES THAT SKEW CORRELATIONS AND
            % EFFECTIVELY FIND SPLITTING EVENTS. NOTE: IN THIS MANNER, THE 'LOW' TRACK IS
            % CONSIDERED THE INITIAL SINGLE TRACK
            
            avg.green.high = getTrackUnique(avg.green.low, avg.green.high_);
            avg.red.high = getTrackUnique(avg.red.low, avg.red.high_);
            
            % DETERMINE SPLIT TIME FROM POINT WHERE HIGH TRACK BEGINS
            
            split.g = getTrackSplit(avg.green.low,avg.green.high);
            split.r = getTrackSplit(avg.red.low, avg.red.high);
            
            % INITIALIZE THE SEGMENTED TRAJECTORIES
            
            g_init = avg.green.low;
            r_init = avg.red.low;
            
            g_low = avg.green.low;
            r_low = avg.red.low;
            
            ss = numel(g_init);
            
            % SPLIT INTO DIFFERNT TRACKS BUT MAINTAIN THE ORIGINAL MATRIX SIZE
            
            for ii = split.g:ss
                g_init(ii) = NaN;
            end
            
            for ii = 1:split.g
                g_low(ii) = NaN;
            end
            
            for ii = split.r:ss
                r_init(ii) = NaN;
            end
            
            for ii = 1:split.r
                r_low(ii) = NaN;
            end
            
            g_high = avg.green.high;
            
            r_high = avg.red.high;

            % GET RELATIVE TRACK COORDS
            
            g_init_rel = getTrackRelative(g_init,L,L0);
            r_init_rel = getTrackRelative(r_init,L,L0);
            g_high_rel = getTrackRelative(g_high,L,L0);
            g_low_rel  = getTrackRelative(g_low,L,L0);
            r_high_rel = getTrackRelative(r_high,L,L0);
            r_low_rel  = getTrackRelative(r_low,L,L0);
            
            
        end

    end
    
    rel.g_init = g_init_rel;
    rel.r_init = r_init_rel;
    rel.g_high = g_high_rel;
    rel.g_low = g_low_rel;
    rel.r_high = r_high_rel;
    rel.r_low = r_low_rel;
    
end


function [track] = getTrackRelative(track,L,L0)

global CONST
PixelSize    = CONST.getLocusTracks.PixelSize;

L_ = L0./L;

right = L/2 - track; 

left = L/2 + track;

track.right = right.*L_;
track.left  = left.*L_;


end


function [L,L0] = getCellLength( data )

global CONST
PixelSize    = CONST.getLocusTracks.PixelSize;

ss = numel( data.CellA );

L0 = PixelSize*data.CellA{1}.length(1)*ones(1,ss);

for ii= 1:ss
    
    L(ii) = PixelSize*data.CellA{ii}.length(1);
    
end

end
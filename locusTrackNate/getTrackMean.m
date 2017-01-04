function [trac, track_mean, variance] = getTrackMean( dirname )

% THIS FUNCTION PICKS OUT THE TRAJECTORIES OF AN ENTIRE TRAJECTORY, AND
% FINDS THE MEAN TRAJECTORY (NOT MSD!)
% IT TREATS PICKS OUT THE PRE-SPLIT VALUES, BUT TREATS THE POST-SPLIT AS
% AN ENTIRE TRACK FROM START TO FINISH
% TO CORRECT FOR THE GROWING CELL, IT SUBTRACTS SCALES THE DATA BY THE
% RELATIVE POSITION OF THE STARTING POINT


global CONST

if isempty( CONST )
    if exist('loadConstantsMine','file');
        loadConstantsMine
    else
        loadConstantsMine
    end
end

TimeStep     = CONST.getLocusTracks.TimeStep/60;
PixelSize    = CONST.getLocusTracks.PixelSize;


dir_list = dir([dirname,filesep,'Cell*.mat']);

num_list_ = numel( dir_list );

for jj_list = 1:num_list_
    
    
    filename = [dirname,filesep,dir_list(jj_list).name];
    
    data = load(filename);
    
    [L,L0] = getCellLength( data );
    
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
        
        dx_g_low = NaN;
        dx_g_high = NaN;
        
    else
        
        if isempty(sort.red.low) || isempty(sort.red.high)
            
            disp('Error: missing red track')
            
            dx_r_low = NaN;
            dx_r_high = NaN;
            
        else
            
            % MOST TRACKS ARE THE SAME, SO WILL GET ONE AVERAGE TRACK PER UPPER AND
            % LOWER ROUTE
            
            avg.green.low = getTrackAvg(sort.green.low);
            avg.green.high = getTrackAvg(sort.green.high);
            avg.red.low = getTrackAvg(sort.red.low);
            avg.red.high = getTrackAvg(sort.red.high);
            
            % REMOVE REPETITIVE VALUES IN TRAJECTORIES THAT SKEW CORRELATIONS AND
            % EFFECTIVELY FIND SPLITTING EVENTS. NOTE: IN THIS MANNER, THE 'LOW' TRACK IS
            % CONSIDERED THE INITIAL SINGLE TRACK
            
            % NOTE THAT THIS PROGRAM KEEPS THE INITIAL REGION IN EACH
            % TRACK
            %
            avg.green.high_ = getTrackUnique(avg.green.low, avg.green.high);
            avg.red.high_ = getTrackUnique(avg.red.low, avg.red.high);
            %
            % DETERMINE SPLIT TIME FROM POINT WHERE HIGH TRACK BEGINS
            
            split.g = getTrackSplit(avg.green.low,avg.green.high_);
            split.r = getTrackSplit(avg.red.low, avg.red.high_);
                        
            % DEV2: JUST RESCALE EVERY COORDINATE WITH L(t)
            %
            
            g_low_time = getTimeDepR0( avg.green.low,L0,L );
            g_high_time = getTimeDepR0( avg.green.high,L0,L );
            r_low_time = getTimeDepR0( avg.red.low,L0,L );
            r_high_time = getTimeDepR0( avg.red.high,L0,L );
            
            g_low_scale = getCellScale( avg.green.low,L0,L);
            g_high_scale = getCellScale( avg.green.high,L0,L);
            r_low_scale = getCellScale( avg.red.low,L0,L);
            r_high_scale = getCellScale( avg.red.high,L0,L);
            
            dx_g_low = g_low_scale - g_low_time;
            dx_g_high = g_high_scale  - g_high_time;
            dx_r_low = r_low_scale - r_low_time;
            dx_r_high = r_high_scale - r_high_time;
            
            % ACHTUNG: COMMENT ONE OF THESE OUT!!!!
            %
            
            % OPTION I: KEEP FULL TRACKS
            % SAVE TRACK VECTORS AS ROWS OF BIG MATRIX TO AVERAGE LATER
            % NOTE: THE HIGH AND LOW TRACKS ALSO INCLUDE THE INITIAL REGION
            
%             ss = size(dx_g_low);
%             
%             dx_g_init(jj_list,1:split.g) = (dx_g_low(1:split.g)-dx_g_low(1));
%             dx_r_init(jj_list,1:split.r) = (dx_r_low(1:split.r)-dx_g_low(1));
%             
%             dx_g_low_(jj_list,1:(ss(2))) = (dx_g_low);
%             dx_g_high_(jj_list,1:(ss(2))) = (dx_g_high);
%             
%             dx_r_low_(jj_list,1:(ss(2))) = (dx_r_low);
%             dx_r_high_(jj_list,1:(ss(2))) = (dx_r_high);
            
            
            % OPTION II: SPLIT INTO INTIAL AND POST SPLIT TRACKS
            
            ss = size(dx_g_low);
            
            dx_g_init(jj_list,1:split.g) = dx_g_low(1:split.g);
            dx_r_init(jj_list,1:split.r) = dx_r_low(1:split.r);
            
            dx_g_low_(jj_list,1:(ss(2)-split.g+1)) = dx_g_low(split.g:end) - dx_g_low(split.g);
            dx_g_high_(jj_list,1:(ss(2)-split.g+1)) = dx_g_high(split.g:end)- dx_g_high(split.g);
            
            dx_r_low_(jj_list,1:(ss(2)-split.r+1)) = dx_r_low(split.r:end) - dx_r_low(split.r);
            dx_r_high_(jj_list,1:(ss(2)-split.r+1)) = dx_r_high(split.r:end) - dx_r_high(split.r);
            
            
            
            
            % SPLIT INTO POST-DETERMINED REGIONS
            %
            % ALERT: THIS IS SOMEWHAT SUBJECTIVE AND VERY SPECIFIC TO
            % INDIVIDUAL TRACKS
            %
            % A and B are the beginning and end of Region II
            % D  is the end of Region III
            
            A = 47;
            B = 90;
            C = 150;
            
            AB = A:B;
            BC = B:C;
            
            AB_ = numel(AB);
            BC_ = numel(BC);
            
            try
                dx_g_h_I(jj_list,1:A) = dx_g_high(1:A).^2;
                dx_g_h_II(jj_list,1:AB_) = (dx_g_high(A:B) - dx_g_high(A)).^2;
                dx_g_h_III(jj_list,1:BC_) = (dx_g_high(B:C) - dx_g_high(B)).^2;
                dx_g_h_IV(jj_list,1:(ss(2)-C)) = (dx_g_high(C+1:end)-dx_g_high(C+1)).^2;
                
            catch
                
                disp('error in region fitting');
            
            end
            
        end
    end
end


ssgi = size(dx_g_init);
ssri = size(dx_r_init);
ssgl = size(dx_g_low_);
ssgh = size(dx_g_high_);
ssrl = size(dx_r_low_);
ssrh = size(dx_r_high_);

trac.gi = getNonZero(dx_g_init);
trac.ri = getNonZero(dx_r_init);
trac.gh = getNonZero(dx_g_high_);
trac.gl = getNonZero(dx_g_low_);
trac.rh = getNonZero(dx_r_high_);
trac.rl = getNonZero(dx_r_low_);

[track_mean.gi,variance.gi] = getTrackMean_(dx_g_init);
[track_mean.ri,variance.ri] = getTrackMean_(dx_r_init);
[track_mean.gh,variance.gh] = getTrackMean_(dx_g_high_);
[track_mean.gl,variance.gl] = getTrackMean_(dx_g_low_);
[track_mean.rh,variance.rh] = getTrackMean_(dx_r_high_);
[track_mean.rl,variance.rl] = getTrackMean_(dx_r_low_);


end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,L0] = getCellLength( data )

global CONST
PixelSize    = CONST.getLocusTracks.PixelSize;

ss = numel( data.CellA );

L0 = PixelSize*data.CellA{1}.length(1)*ones(1,ss);

for ii= 1:ss
    
    L(ii) = PixelSize*data.CellA{ii}.length(1);
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [track_] = getTimeDepR0( track,L0,L )

ss = find(~isnan(track),1,'first');

sss = size (track);

r0(1:sss(2)) = track(ss);

% L_ = L./L0;
%
% track_ = r0.*L_;

track_ = r0./L0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [track_] = getCellScale(track,L0,L);

track_ = track./L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


function [msd, variance, err] = getTrackMean_(matrix)

ss = size(matrix);
sss = numel(matrix);

for mm = 1:sss
   
    if matrix(mm) == 0
        matrix(mm) = NaN;
    end
    
end


for ii = 1:ss(2)
    
    msd_ = matrix(:,ii);
    msd__ = msd_(~isnan(msd_));
    
    if numel(msd__) <= 3
        
        err(ii) = NaN;
        msd(ii) = NaN;
        variance(ii) = NaN;
        
    else
        
        err(ii) = std(msd__)/sqrt(numel(msd__));
        msd(ii) = mean(msd__);
        variance(ii) = var(msd__);
        
    end
    
end

% err = err;
% msd = msd;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matrix_] = getNonZero(matrix)

sss = size(matrix);

for ii = 2:sss(2)
    for jj = 1:sss(1)
        if matrix(jj,ii) == 0
            matrix(jj,ii) = NaN;
        end
    end
    
end

matrix_ = matrix;

end

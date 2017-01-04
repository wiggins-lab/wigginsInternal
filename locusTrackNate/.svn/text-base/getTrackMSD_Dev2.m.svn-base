function [msd, variance, err] = getTrackMSD_Dev2( dirname )

%% NJK, MAY 2011
% THIS FUNCTION CALCUATES THE MEAN SQUARED DISPLACEMENT OF THE LOCI. IT
% AUTOMATICALLY SPLITS THE TRACK INTO PRE-SPLIT AND POST-SPLIT REGIONS, AND
% HAS THE OPTION OF EITHER TREATING THE SISTER TRACKS AS STARTING FROM THE
% BEGINNING OF THE TRAJECTORY OR AT THE INSTANT AFTER THE SPLIT.

%% LOAD CONSTANTS
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

%% DETERMINE FILE LIST

dir_list = dir([dirname,filesep,'Cell*.mat']);

num_list_ = numel( dir_list );

%% BEGIN LOOP
for jj_list = 1:num_list_
    
    
    filename = [dirname,filesep,dir_list(jj_list).name];
    
    data = load(filename);

    %% GET CELL LENGTH AS A FUNCTION OF TIME
    [L,L0] = getCellLength( data );
    
    %% IF TRACKS HAVE ALREADY BEEN SAVED, USE THOSE
    if ~isfield(data,'track1') || ~isfield(data,'track2')
        
        disp([filename, ' : Getting New Locus Tracks']);
        tracks = getLocusTracksDev_lite(data);
        
    else
        
        disp([filename, ' : Using Existing Locus Tracks']);
        tracks = data;
        
    end
    

    %% SORT TRACKS IN XTRACK INTO UPPER AND LOWER TRACKS, DETERMINED BY THE
    % SIGN OF THE LAST 5 POSITION COORDINATES (NEGATIVE ==> LOWER)
    
    sort.green  = getTrackSort(tracks.track1.xtrack.*PixelSize);
    sort.red = getTrackSort(tracks.track2.xtrack.*PixelSize);
    
    %% MAKE SURE THAT THERE ARE TWO OF EACH TRACK
    
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
            
            %% MOST TRACKS ARE THE SAME, SO WILL GET ONE AVERAGE TRACK PER UPPER AND
            % LOWER ROUTE
            
            avg.green.low = getTrackAvg(sort.green.low);
            avg.green.high = getTrackAvg(sort.green.high);
            avg.red.low = getTrackAvg(sort.red.low);
            avg.red.high = getTrackAvg(sort.red.high);
            
            %% REMOVE REPETITIVE VALUES IN TRAJECTORIES THAT SKEW CORRELATIONS AND
            % EFFECTIVELY FIND SPLITTING EVENTS. NOTE: IN THIS MANNER, THE 'LOW' TRACK IS
            % CONSIDERED THE INITIAL SINGLE TRACK
            %
            % NOTE: THIS IS ONLY USED TO FIND THE SPLIT POINT IN OPTION I
            % BELOW
            
            avg.green.high_ = getTrackUnique(avg.green.low, avg.green.high);
            avg.red.high_ = getTrackUnique(avg.red.low, avg.red.high);
            
            
            %% DETERMINE SPLIT TIME FROM POINT WHERE HIGH TRACK BEGINS
            
            split.g = getTrackSplit(avg.green.low,avg.green.high_);
            split.r = getTrackSplit(avg.red.low, avg.red.high_);
                       
            %% RESCALE EVERY COORDINATE WITH L(t) TO FACTOR OUT CELL GROWTH
            % IN TRAJECTORY
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

            %% MAKE TRACKS
            %
%           % ALERT: COMMENT ONE OF THESE OPTIONS OUT OR BE DOOMED            
%           
%             
%             % OPTION I: KEEP FULL LENGTH TRACK
%             % SAVE TRACK VECTORS AS ROWS OF BIG MATRIX TO AVERAGE LATER
%             % NOTE: THE HIGH AND LOW TRACKS ALSO INCLUDE THE INITIAL REGION
%             %
%             
%             
%             ss = size(dx_g_low);
%             
%             dx_g_init(jj_list,1:split.g) = (dx_g_low(1:split.g)-dx_g_low(1)).^2;
%             dx_r_init(jj_list,1:split.r) = (dx_r_low(1:split.r)-dx_g_low(1)).^2;
%             
%             dx_g_low_(jj_list,1:(ss(2))) = (dx_g_low).^2;
%             dx_g_high_(jj_list,1:(ss(2))) = (dx_g_high).^2;
%             
%             dx_r_low_(jj_list,1:(ss(2))) = (dx_r_low).^2;
%             dx_r_high_(jj_list,1:(ss(2))) = (dx_r_high).^2;
            
            % OPTION II: SPLIT INTO INTIAL AND POST SPLIT TRACKS
            
            ss = size(dx_g_low);
            
            dx_g_init(jj_list,1:split.g) = dx_g_low(1:split.g).^2;
            dx_r_init(jj_list,1:split.r) = dx_r_low(1:split.r).^2;
            
            dx_g_low_(jj_list,1:(ss(2)-split.g+1)) = (dx_g_low(split.g:end) - dx_g_low(split.g)).^2;
            dx_g_high_(jj_list,1:(ss(2)-split.g+1)) = (dx_g_high(split.g:end)- dx_g_high(split.g)).^2;
            
            dx_r_low_(jj_list,1:(ss(2)-split.r+1)) = (dx_r_low(split.r:end) - dx_r_low(split.r)).^2;
            dx_r_high_(jj_list,1:(ss(2)-split.r+1)) = (dx_r_high(split.r:end) - dx_r_high(split.r)).^2;
            
            
            %% OPTIONAL: SPLIT INTO POST-DETERMINED REGIONS
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

%% FIND SIZE OF TRAJECTORY MATRICES
ssgi = size(dx_g_init);
ssri = size(dx_r_init);
ssgl = size(dx_g_low_);
ssgh = size(dx_g_high_);
ssrl = size(dx_r_low_);
ssrh = size(dx_r_high_);

%% DETERMINE MEAN AND VARIANCE OF SQUARED DISPLACEMENT 
[msd_gi, var_gi, err_gi] = getMSDavg(dx_g_init);
[msd_ri, var_ri, err_ri] = getMSDavg(dx_r_init);
[msd_gl, var_gl, err_gl] = getMSDavg(dx_g_low_);
[msd_gh, var_gh, err_gh] = getMSDavg(dx_g_high_);
[msd_rl, var_rl, err_rl] = getMSDavg(dx_r_low_);
[msd_rh, var_rh, err_rh] = getMSDavg(dx_r_high_);

[msd_ghI, err_ghI] = getMSDavg(dx_g_h_I);
[msd_ghII, err_ghII] = getMSDavg(dx_g_h_II);
[msd_ghIII, err_ghIII] = getMSDavg(dx_g_h_III);
[msd_ghIV, err_ghIV] = getMSDavg(dx_g_h_IV);

msd.gi = msd_gi;
msd.ri = msd_ri;
msd.gl = msd_gl;
msd.gh = msd_gh;
msd.rl = msd_rl;
msd.rh = msd_rh;

msd.gI = msd_ghI;
msd.gII = msd_ghII;
msd.gIII = msd_ghIII;
msd.gIV = msd_ghIV;

variance.gi = var_gi;
variance.ri = var_ri;
variance.gl = var_gl;
variance.gh = var_gh;
variance.rl = var_rl;
variance.rh = var_rh;

err.gi = err_gi;
err.ri = err_ri;
err.gl = err_gl;
err.gh = err_gh;
err.rl = err_rl;
err.rh = err_rh;
err.gI = err_ghI;
err.gII = err_ghII;
err.gIII = err_ghIII;
err.gIV = err_ghIV;

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,L0] = getCellLength( data )

%% FUNCTION TO GET CELL LENGTH AS A FUNTION OF TIME

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

%% FUNCTION TO DETERMINE THE TIME DEPENDENT `INITIAL' COORDINATE OF EACH
% TRJECTORY
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

%% FUNCTION TO SCALE THE TRAJECTROY BY CELL LENGTH

track_ = track./L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


function [msd, variance, err] = getMSDavg(matrix)

%% FUNCTION TO GET MEAN AND VARIANCE FROM SQUARED TRAJECTORY MATRICES
ss = size(matrix);

for ii = 1:ss(2)
    
    msd_ = matrix(:,ii);
    msd__ = msd_(msd_ > 0);
    
    if numel(msd__) <= 3
        
        err(ii) = NaN;
        msd(ii) = NaN;
        
    else
        
        err(ii) = std(msd__)/sqrt(numel(msd__));
        msd(ii) = mean(msd__);
        variance(ii) = var(msd__);
        
    end
    
end

% err = err;
% msd = msd;

end

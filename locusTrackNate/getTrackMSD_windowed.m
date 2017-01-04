function [slope] = getTrackMSD_windowed( dirname, window )

% THIS FUNCTION IS DESIGNED TO ULTIMATELY PRODUCE A PLOT
% OF THE FIT PARAMETER alpha ON A LOGLOG PLOT AS A FUNCTION 
% OF TIME FOR A PRESET WINDOW OF POINTS.  IT FIRST WINDOWS THE TRACK
% STRUCTURE, THEN COMPUTES MSD FOR THAT WINDOW (STARTING FROM ZERO) 
% THEN FITS THE LOGLOG PLOT OF THAT WINDOWED MSD.



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


dir_list = dir([dirname,filesep,'Cell*.mat']);

num_list_ = numel( dir_list );

dx_g_low_ = [];
dx_g_high_ = [];
dx_r_low_ = [];
dx_r_high_ = [];


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
            
            %avg.green.high = getTrackUnique(avg.green.low, avg.green.high_);
            %avg.red.high = getTrackUnique(avg.red.low, avg.red.high_);
            
            % DETERMINE SPLIT TIME FROM POINT WHERE HIGH TRACK BEGINS
            
            split.g = getTrackSplit(avg.green.low,avg.green.high);
            split.r = getTrackSplit(avg.red.low, avg.red.high);

            ss = numel(avg.green.low);
            
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
            
            dx_g_low(jj_list,1:ss) = g_low_scale - g_low_time;
            dx_g_high(jj_list,1:ss) = g_high_scale  - g_high_time;
            dx_r_low(jj_list,1:ss) = r_low_scale - r_low_time;
            dx_r_high(jj_list,1:ss) = r_high_scale - r_high_time;
            
        end
    end
end

% SIZE OF DISPLACEMENT MATRICES

ssgl = size(dx_g_low);
ssgh = size(dx_g_high);
ssrl = size(dx_r_low);
ssrh = size(dx_r_high);

slope.gh = getMSDslope(dx_g_high,window);            
slope.gl = getMSDslope(dx_g_low,window);            
slope.rh = getMSDslope(dx_r_high,window);            
slope.rl = getMSDslope(dx_r_low,window);


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

% % FOR DEV
% 
% L_ = L./L0;
% 
% track_ = r0.*L_;

% FOR DEV2

track_ = r0./L0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [track_] = getCellScale(track,L0,L);

track_ = track./L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function slope = getMSDslope( matrix, window );

% THIS FUNCTION WINDOWS THE TRACK STRUCTURE TO CREATE SHORT MSD SEGMENTS
% THAT ALL START FROM ZERO FOR LOGLOG MSD ANALYSIS
%
%  

ss = size(matrix);
win = floor(window/2);

len = floor(ss(2)/win);

for ii = 1:(len-2)
    
    iii = ii*win + 1;
    kk = iii-win;
    mm = kk + window;
    
%for ii = 1:(len-1)
    
%     kk = ii*window - window +1;
%     mm = kk + window;
        
    mat = matrix(:,(kk:mm));
    
    for jj = 1:window    
        mat_(:,jj) = mat(:,jj) - mat(:,1);
    end
    
    msd = getMSDavg(mat_);
    
    A = find(~isnan(msd),1,'first');
    B = find(~isnan(msd),1,'last');
    
% %     % METHOD 1: FINITE DIFFERENCES
% %     
%     rise = log10(msd(B))- log10(msd(A));
%     run = log10(B-A);
%     
%     slope(iii) = rise/run;
%     
%     % METHOD 2: LINEAR FIT USING polyfit
%     
    ln_msd = log10(msd);
    
    t = 1:numel(msd);
    
    ln_t = log10(t);
    
    [P,s] = polyfit(ln_t(~isnan(ln_msd)),ln_msd(~isnan(ln_msd)),1);
    
    slope(ii) = P(1);
        
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function msd = getMSDavg(matrix)

ss = size(matrix);

for ii = 1:ss(2)
    
    msd_ = matrix(:,ii);
    msd__ = msd_(msd_ > 0);
    
    if numel(msd__) <= 3

        msd(ii) = NaN;

    else
        
    msd(ii) = mean(msd__);
        
    end
    
end

end

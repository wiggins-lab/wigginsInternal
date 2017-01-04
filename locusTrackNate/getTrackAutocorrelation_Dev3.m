function [tots] = getTrackAutocorrelation_Dev3( filename, len, view_flag )

%% THIS FUNCTION IS USED TO FIND THE AUTOCORRELATION (VELOCITY-VELOCITY
% CORRELATIONS) IN A SINGLE TRAJECTORY.  IT IS CALLED FROM
% BatchgetTrackAuto.m TO PROCESS AN ENTIRE TRAJECORY
%
% NJK

%% INITIALIZE CONSTANTS
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

%% DETERMINE THE LENGTH OF CORRELATED MOTION IN QUESTION
% LEN IS THE LENGTH OF AUTOCORRELATION, 2 IS NEXT TIME STEP

if nargin < 2 || isempty( len )
    
    len = 2;
    
end

% DEFAULT BEHAVIOR IS SHOW SUBPLOT
if nargin < 3 || isempty ( view_flag )
    
    view_flag = 0;
    
end

%% LOAD CELL DATA

data = load(filename);

% GET TRACK STRUCTURE, INCLUDING XTRACK AND SPLITS
if ~isfield(data,'track1') || ~isfield(data,'track2')
    
    disp([filename, ' : Getting New Locus Tracks']);
    tracks = getLocusTracksDev_lite(data);
    
else
    
    disp([filename, ' : Using Existing Locus Tracks']);
    tracks = data;
    
end

sort.green = [];
sort.red = [];

%% SORT TRACKS IN XTRACK INTO UPPER AND LOWER TRACKS, DETERMINED BY THE
% SIGN OF THE LAST 5 POSITION COORDINATES (NEGATIVE ==> LOWER)

sort.green  = getTrackSort(tracks.track1.xtrack.*PixelSize);
sort.red = getTrackSort(tracks.track2.xtrack.*PixelSize);

%% MAKE SURE THAT THERE ARE TWO OF EACH TRACK

if isempty(sort.green.low) || isempty(sort.green.high)
    
    disp('Error: missing green track')
    
    auto_g_init = NaN;
    auto_r_init = NaN;
    
    auto_g_low = NaN;
    auto_g_high = NaN;
    auto_r_low = NaN;
    auto_r_high = NaN;
    
else
    
    if isempty(sort.red.low) || isempty(sort.red.high)
        
        disp('Error: missing red track')
        
        auto_g_init = NaN;
        auto_r_init = NaN;
        
        auto_g_low = NaN;
        auto_g_high = NaN;
        auto_r_low = NaN;
        auto_r_high = NaN;
        
        
    else
        
        
        %% MOST TRACKS ARE THE SAME, SO WILL GET ONE AVERAGE TRACK PER UPPER AND
        % LOWER ROUTE
        
        avg.green.low = getTrackAvg(sort.green.low);
        avg.green.high_ = getTrackAvg(sort.green.high);
        avg.red.low = getTrackAvg(sort.red.low);
        avg.red.high_ = getTrackAvg(sort.red.high);
        
        %% REMOVE REPETITIVE VALUES IN TRAJECTORIES THAT SKEW CORRELATIONS AND
        % EFFECTIVELY FIND SPLITTING EVENTS. NOTE: IN THIS MANNER, THE 'LOW' TRACK IS
        % CONSIDERED THE INITIAL SINGLE TRACK
        
        avg.green.high = getTrackUnique(avg.green.low, avg.green.high_);
        avg.red.high = getTrackUnique(avg.red.low, avg.red.high_);
        
        %% DETERMINE SPLIT TIME FROM POINT WHERE HIGH TRACK BEGINS
        
        split.g = getTrackSplit(avg.green.low,avg.green.high);
        split.r = getTrackSplit(avg.red.low, avg.red.high);
        
        
        
        %% GET DISPLACMENTS BETWEEN SUCCESIVE TIME POINTS
        
        dx_g_low = getTrackDisplacements(avg.green.low') ;
        dx_g_high = getTrackDisplacements(avg.green.high');
        
        dx_r_low = getTrackDisplacements(avg.red.low') ;
        dx_r_high = getTrackDisplacements(avg.red.high');
        
        
        %% GET AUTOCORRELATION (FUNCTION AT BOTTOM OF PROGRAM)
        % THIS OUTPUTS A FORWARD AND BACKWARD CALCULATION
        
        auto_g_init = getAutoDev(dx_g_low(1:split.g));
        auto_r_init = getAutoDev(dx_r_low(1:split.r));
        
        auto_g_low = getAutoDev(dx_g_low(split.g:end));
        auto_g_high = getAutoDev(dx_g_high);
        auto_r_low = getAutoDev(dx_r_low(split.r:end));
        auto_r_high = getAutoDev(dx_r_high);
        
        %% FLIP AROUND THE BACKWARD TRACE FOR PLOTTING
        
        [x_g_ib, x_g_if] = getX(auto_g_init);
        [x_r_ib, x_r_if] = getX(auto_r_init);
        [x_g_lb, x_g_lf] = getX(auto_g_low);
        [x_g_hb, x_g_hf] = getX(auto_g_high);
        [x_r_lb, x_r_lf] = getX(auto_r_low);
        [x_r_hb, x_r_hf] = getX(auto_r_high);
        
        
        if view_flag == 0
            
            figure(10);
            hold on          
            
            subplot(3,2,1)
            plot(x_g_if, auto_g_init.forward)
            hold on
            plot(x_g_ib, auto_g_init.back);
            title('g init');
                        
            subplot(3,2,2)
            plot(x_r_if, auto_r_init.forward)
            hold on
            plot(x_r_ib, auto_r_init.back);
            title('r init');
            
            subplot(3,2,3)
            plot(x_g_lf, auto_g_low.forward)
            hold on
            plot(x_g_lb, auto_g_low.back);
            title('g low');
            
            subplot(3,2,4)
            plot(x_g_hf, auto_g_high.forward)
            hold on
            plot(x_g_hb, auto_g_high.back);
            title('g high');
                        
            subplot(3,2,5)
            plot(x_r_lf, auto_r_low.forward)
            hold on
            plot(x_r_lb, auto_r_low.back);
            title('r low');
            
            subplot(3,2,6)
            plot(x_r_hf, auto_r_high.forward)
            hold on
            plot(x_r_hb, auto_r_high.back);
            title('r high');
            
        end
    end
end

%% EXPORT MATRICES TO BatchGetTrackAuto.m

tots.auto_g_init_tot = auto_g_init;
tots.auto_r_init_tot = auto_r_init;

tots.auto_g_low_tot = auto_g_low;
tots.auto_g_high_tot =  auto_g_high;
tots.auto_r_low_tot =  auto_r_low;
tots.auto_r_high_tot = auto_r_high;

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Auto] = getAutoDev(track)

%% THIS FUNCTION GET AUTOCORRELTAITON VALUES 'AHEAD' AND BEHIND THE POINT
% IN QUESTION, CONVENIENT FOR PLOTTING
ss = size(track);

smid_tot= ceil(ss(1)/2);

for jj = 2:smid_tot
    
    for ii = 0:(jj-1)
        
        Auto_for(jj,ii+1) = track(jj)*track(jj+ii);
        
        Auto_back(jj,ii+1) = track(jj)*track(jj-ii);
        
    end
    
end

for jj = smid_tot:(ss(1)-1)
    
    for ii = 0:(ss(1)-jj-1)
        
        Auto_for(jj,ii+1) = track(jj)*track(jj+ii);
        
        Auto_back(jj,ii+1) = track(jj)*track(jj-ii);
        
    end
    
end


ss_f = size(Auto_for);
ss_b = size(Auto_back);

for kk = 1:ss_f(2)
    
    Auto_for_avg(kk) = mean(Auto_for(~isnan(Auto_for(:,kk)),kk));
    
end

for kk = 1:ss_b(2)
    
    Auto_back_avg(kk) = mean(Auto_back(~isnan(Auto_back(:,kk)),kk));
    
end

Auto.forward = Auto_for_avg;
Auto.back = Auto_back_avg;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xb, xf] = getX(track)

%% USED FOR PLOTTING (I THINK)
lenb = numel(track.back)-1;

xb = 0:-1:-lenb;

lenf = numel(track.forward)-1;

xf = 0:1:lenf;

end


% OLD CORR CALCULATION
%
% function auto = getTrackAuto(track,len)
% 
% ssA = size(track);
% 
% auto = track(len:end,:).*track(1:(ssA(1)-len+1),:);
% 
% end
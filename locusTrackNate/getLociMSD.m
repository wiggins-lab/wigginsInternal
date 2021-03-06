function [msd, err] = getLociMSD( dirname )

% THIS FUNCTION CALCUATES THE MEAN SQUARED DISPLACEMENT BETWEEN THE LOCI,
% IT PICKS OUT THE PRE-SPLIT VALUES, BUT TREATS THE POST-SPLIT AS
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

% figure(10);
% clf;
% figure(10);
% hold on;
%
% figure(11);
% clf;
% figure(11);
% hold on;
%
% figure(12);
% clf
% figure(12);
% hold on


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
            
            % NOTE THAT THIS PROGRAM KEEPS THE INITIAL REGION IN EACH
            % TRACK
            %
            avg.green.high_ = getTrackUnique(avg.green.low, avg.green.high);
            avg.red.high_ = getTrackUnique(avg.red.low, avg.red.high);
            %
            % DETERMINE SPLIT TIME FROM POINT WHERE HIGH TRACK BEGINS
            
            split.g = getTrackSplit(avg.green.low,avg.green.high_);
            split.r = getTrackSplit(avg.red.low, avg.red.high_);
            
            % DEV: INSEAD OF THE INITIAL POSITION STAYING CONSTANT,
            % WE MAY NEED TO SUBTRACT AN ADJUSTED DISTANCE DUE TO INFLATION
            %
            % g_low_time = getTimeDepR0( avg.green.low,L0,L );
            % g_high_time = getTimeDepR0( avg.green.high,L0,L );
            % r_low_time = getTimeDepR0( avg.red.low,L0,L );
            % r_high_time = getTimeDepR0( avg.red.high,L0,L );
            %
            %
            % dx_g_low = avg.green.low - g_low_time;
            % dx_g_high = avg.green.high - g_high_time;
            % dx_r_low = avg.red.low - r_low_time;
            % dx_r_high = avg.red.high - r_high_time;
            
            
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
            
            % SAVE TRACK VECTORS AS ROWS OF BIG MATRIX TO AVERAGE LATER
            % NOTE: THE HIGH AND LOW TRACKS ALSO INCLUDE THE INITIAL REGION
            
            sg = size(dx_g_low);
            sr = size(dx_r_low);                        
            
            sg_ = sg(2) - split.g +1;
            sr_ = sr(2) - split.r +1;
           
            sd_g12(jj_list,1:sg_) = (dx_g_high(split.g:end) - dx_g_low(split.g:end)).^2;
            sd_r12(jj_list,1:sr_)  = (dx_r_high(split.r:end) - dx_r_low(split.r:end)).^2;
            
            
%             % SPLIT INTO POST-DETERMINED REGIONS
%             %
%             % ALERT: THIS IS SOMEWHAT SUBJECTIVE AND VERY SPECIFIC TO
%             % INDIVIDUAL TRACKS
%             %
%             % A and B are the beginning and end of Region II
%             % D  is the end of Region III
%             
%             A = 47;
%             B = 90;
%             C = 150;
%             
%             AB = A:B;
%             BC = B:C;
%             
%             AB_ = numel(AB);
%             BC_ = numel(BC);
%             
%             try
%                 dx_g_h_I(jj_list,1:A) = dx_g_high(1:A).^2;
%                 dx_g_h_II(jj_list,1:AB_) = (dx_g_high(A:B) - dx_g_high(A)).^2;
%                 dx_g_h_III(jj_list,1:BC_) = (dx_g_high(B:C) - dx_g_high(B)).^2;
%                 dx_g_h_IV(jj_list,1:(ss(2)-C)) = (dx_g_high(C+1:end)-dx_g_high(C+1)).^2;
%                 
%             catch
%                 
%                 disp('error in region fitting');
%             
%             end
            
        end
    end
end

[msd_g12, err_g12] = getMSDavg(sd_g12);
[msd_r12, err_r12] = getMSDavg(sd_r12);

msd.g12 = msd_g12;
msd.r12 = msd_r12;

err.g12 = err_g12;
err.r12 = err_r12;


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


function [msd, err] = getMSDavg(matrix)

ss =  size(matrix);

for ii = 1:ss(2)
    
    msd_ = matrix(:,ii);
    msd__ = msd_(msd_ > 0);
    
    if numel(msd__) <= 3
        
        err(ii) = NaN;
        msd(ii) = NaN;
        
    else
        
        err(ii) = std(msd__)/sqrt(numel(msd__));
        msd(ii) = mean(msd__);
        
    end
    
end

% err = err;
% msd = msd;

end
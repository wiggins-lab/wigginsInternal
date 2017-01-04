function [corr_sum] = getTrackCorrelation_Dev2( filename, jj )

% NJK, MARCH 2011.  THIS FUNCTION IS USED TO DETERMINE CORRELATIONS BETWEEN
% THE MOVEMENT OF THE LOCI.  IT CALCUATES DISPLACMENT VECTORS AT EACH TIME
% POINT FOR EACH LOCI, THEN USES A NORMALIZED VECTOR PRODUCT BETWEEN
% DISPLACMENT VECTORS AS THE MEASURE OF CORRELATION.

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


data = load(filename);

% GET TRACK STRUCTURE, INCLUDING XTRACK AND SPLITS

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
    
    corr_sum.init = NaN;
    corr_sum.g1g2 = NaN;
    corr_sum.r1r2 = NaN;
    corr_sum.g1r1 = NaN;
    corr_sum.g1r2 = NaN;
    corr_sum.g2r1 = NaN;
    corr_sum.g2r2 = NaN;
    
else
    
if isempty(sort.red.low) || isempty(sort.red.high)
    
    disp('Error: missing red track')
    
    corr_sum.init = NaN;
    corr_sum.g1g2 = NaN;
    corr_sum.r1r2 = NaN;
    corr_sum.g1r1 = NaN;
    corr_sum.g1r2 = NaN;
    corr_sum.g2r1 = NaN;
    corr_sum.g2r2 = NaN;

else
    

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


% GET DISPLACMENTS BETWEEN SUCCESIVE TIME POINTS

dx_g_low = getTrackDisplacements(avg.green.low') ;
dx_g_high = getTrackDisplacements(avg.green.high');

dx_r_low = getTrackDisplacements(avg.red.low') ;
dx_r_high = getTrackDisplacements(avg.red.high');

% CALCULATE THE MAGNITUDE OF EACH TRAJECTORY FOR NORALIZATION

mag_g_low = sqrt(sum(dx_g_low(~isnan(dx_g_low)).^2));
mag_g_high = sqrt(sum(dx_g_high(~isnan(dx_g_high)).^2));
mag_r_low = sqrt(sum(dx_r_low(~isnan(dx_r_low)).^2));
mag_r_high = sqrt(sum(dx_r_high(~isnan(dx_r_high)).^2));

% CALCULATE CORRELATIONS. EFFECTIVELY DOT PRODUCTS (WITHOUT SUMMING)
% BETWEEN DISPLACEMENT VECTORS AT EACH TIME POINT, NORMALIZED TO THE
% LENGTH OF THE WHOLE TRAJECTORY

corr_g1g2 = dx_g_low.*dx_g_high;
corr_r1r2 = dx_r_low.*dx_r_high;
corr_g1r1 = dx_g_low.*dx_r_low;
corr_g1r2 = dx_g_low.*dx_r_high;
corr_g2r1 = dx_g_high.*dx_r_low;
corr_g2r2 = dx_g_high.*dx_r_high;

% SUM UP THE CORRELATION VECTORS OF THE INITIAL TRACKS AND POST-SPIIT
% TRACKS

init = corr_g1r1(1:split.r);
post = corr_g1r1(split.g:end);

corr_sum.init = sum(init(~isnan(init)));
corr_sum.g1g2 = sum(corr_g1g2(~isnan(corr_g1g2)));
corr_sum.r1r2 = sum(corr_r1r2(~isnan(corr_r1r2)));
corr_sum.g1r1 = sum(post(~isnan(post)));
corr_sum.g1r2 = sum(corr_g1r2(~isnan(corr_g1r2)));
corr_sum.g2r1 = sum(corr_g2r1(~isnan(corr_g2r1)));
corr_sum.g2r2 = sum(corr_g2r2(~isnan(corr_g2r2)));

% PLOT HISTOGRAMS OF CORRELATION MAGNITUDES 

figure(10);
hold on


subplot(4,2,1)
hist(corr_g1r1(1:split.g))
title('G1-R1 initial')
hold on

subplot(4,2,3)
hist(corr_g1r1(split.g:end))
title('G1-R1 post-split')
hold on

subplot(4,2,4)
hist(corr_g2r2)
title('G2-R2')
hold on

subplot(4,2,5)
hist(corr_g1g2)
title('G1-G2')
hold on

subplot(4,2,6)
hist(corr_r1r2)
title('R1-R2')
hold on

subplot(4,2,7)
hist(corr_g1r2)
title('G1-R2')
hold on

subplot(4,2,8)
hist(corr_g2r1)
title('G2-R1')
hold on


    
end
end

end









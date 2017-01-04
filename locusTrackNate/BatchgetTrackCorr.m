function [] = BatchgetTrackCorr( dirname )

%% THIS FUNCTION IS USED TO GET CORRELATION FUNCTIONS FOR AN ENTIRE
% DIRECTORY.  IT USES getTrackCorrelation_Dev2.m.
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

figure(10);
clf
figure(11);
clf

%% READ FILES
dir_list = dir([dirname,filesep,'Cell*.mat']);

num_list_ = numel( dir_list );

corr_sum.init = [];
corr_sum.g1r1 = [];
corr_sum.g1g2 = [];
corr_sum.r1r2 = [];
corr_sum.g1r2 = [];
corr_sum.g2r1 = [];
corr_sum.g2r2 = [];

%% GET CORRELATIONS
for jj = 1:num_list_
    
    filename = [dirname,filesep,dir_list(jj).name];
    
    corr_sum_ = getTrackCorrelation_Dev2(filename,jj);
       
    corrs.init = corr_sum_.init;
    corrs.g1r1 = corr_sum_.g1r1;
    corrs.g1g2 = corr_sum_.g1g2;
    corrs.r1r2 = corr_sum_.r1r2;
    corrs.g1r2 = corr_sum_.g1r2;
    corrs.g2r1 = corr_sum_.g2r1;
    corrs.g2r2 = corr_sum_.g2r2;
    
    corr_sum.init(jj) = corr_sum_.init;
    corr_sum.g1r1(jj) = corr_sum_.g1r1;
    corr_sum.g1g2(jj) = corr_sum_.g1g2;
    corr_sum.r1r2(jj) = corr_sum_.r1r2;
    corr_sum.g1r2(jj) = corr_sum_.g1r2;
    corr_sum.g2r1(jj) = corr_sum_.g2r1;
    corr_sum.g2r2(jj) = corr_sum_.g2r2;
    
    
end


%% FIND THE SUM OF CORRELATIONS (TO DETERMINE POSITIVE, NEGATIVE, OR NO
% CORRELATIONS)
init_tot = sum(corr_sum.init(~isnan(corr_sum.init)));
g1r1_tot = sum(corr_sum.g1r1(~isnan(corr_sum.g1r1)));
g1g2_tot = sum(corr_sum.g1g2(~isnan(corr_sum.g1g2)));
r1r2_tot = sum(corr_sum.r1r2(~isnan(corr_sum.r1r2)));
g1r2_tot = sum(corr_sum.g1r2(~isnan(corr_sum.g1r2)));
g2r1_tot = sum(corr_sum.g2r1(~isnan(corr_sum.g2r1)));
g2r2_tot = sum(corr_sum.g2r2(~isnan(corr_sum.g2r2)));

%% PLOT RESULTS AS HISTOGRAMS

figure(11);

subplot(4,2,1)
hist(corr_sum.init);
title(['Initial = ', num2str(init_tot)])

subplot(4,2,3)
hist(corr_sum.g1r1)
title(['G1-R1 = ', num2str(g1r1_tot)])

subplot(4,2,4)
hist(corr_sum.g2r2)
title(['G2-R2 = ', num2str(g2r2_tot)])

subplot(4,2,5)
hist(corr_sum.g1g2)
title(['G1-G2 = ', num2str(g1g2_tot)])

subplot(4,2,6)
hist(corr_sum.r1r2)
title(['R1-R2 = ', num2str(r1r2_tot)])

subplot(4,2,7)
hist(corr_sum.g1r2)
title(['G1-R2 = ', num2str(g1r2_tot)])

subplot(4,2,8)
hist(corr_sum.g2r1)
title(['G2-R1 = ', num2str(g2r1_tot)])





end
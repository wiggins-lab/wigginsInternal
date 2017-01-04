function [avgs_corr] = BatchgetTrackAuto( dirname, corr_len, view_flag )

%% THIS FUNCTION IS USED TO GET AUTOCORRELATION GRAPHS FOR AN ENTIRE
% TRAJECTORY, USING THE PROGRAM getTrackAutocorrelation_Dev3.m

%% INITIALIZE CONSTANTS
global CONST

if isempty( CONST )
    if exist('loadConstantsMine','file');
        loadConstantsMine
    else
        loadConstants
    end
end

%% DETERMINE LENGTH OF CORRELATIONS
if nargin < 2 || isempty( corr_len )
    
    corr_len = 1;

end

%% SURPRESS OUTPUT IN getTrackAutocorrelation_Dev3.m BY DEFAULT
if nargin < 3 || isempty ( view_flag )
    
    view_flag = 0;
    
end


TimeStep     = CONST.getLocusTracks.TimeStep/60;
PixelSize    = CONST.getLocusTracks.PixelSize;

figure(10);
clf
figure(16);
clf

dir_list = dir([dirname,filesep,'Cell*.mat']);

num_list_ = numel( dir_list );

corr_length = corr_len + 1;

%% BEGIN LOOP TO GET AUTOCORRELATIONS
for jj = 1:num_list_
    
    filename = [dirname,filesep,dir_list(jj).name];
    
    tots = getTrackAutocorrelation_Dev3( filename, corr_length, view_flag);
    
    g_init{jj} = tots.auto_g_init_tot;
    r_init{jj} = tots.auto_r_init_tot;
    g_low{jj} = tots.auto_g_low_tot;
    g_high{jj} = tots.auto_g_high_tot;
    r_low{jj} = tots.auto_r_low_tot;
    r_high{jj} = tots.auto_r_high_tot;

end

%% GET AVERAGE VALUES
    avg_g_init = getAutoAvg(g_init);
    avg_r_init = getAutoAvg(r_init);
    avg_g_low = getAutoAvg(g_low);
    avg_g_high = getAutoAvg(g_high);
    avg_r_low = getAutoAvg(r_low);
    avg_r_high = getAutoAvg(r_high);
    
    avgs_corr.g_init = avg_g_init;
    avgs_corr.r_init = avg_r_init;
    avgs_corr.g_low = avg_g_low;
    avgs_corr.g_high = avg_g_high;
    avgs_corr.r_low = avg_r_low;
    avgs_corr.r_high = avg_r_high;
    
    sum_g_init = corr_sums(avg_g_init);
    sum_r_init = corr_sums(avg_r_init);
    sum_g_low = corr_sums(avg_g_low);
    sum_g_high = corr_sums(avg_g_high);
    sum_r_low = corr_sums(avg_r_low);
    sum_r_high = corr_sums(avg_r_high);

    %% PLOT RESULTS
    
    figure(16);
    hold on
    
    subplot(3,2,1)
    plot(0:15, avg_g_init.forward(1:16),'-og')
    hold on
    plot(0:-1:-15, avg_g_init.back(1:16),'-og');
    title(['g init =', num2str(sum_g_init)]);
    
    subplot(3,2,2)
    plot(0:15, avg_r_init.forward(1:16),'-or')
    hold on
    plot(0:-1:-15, avg_r_init.back(1:16),'-or');
    title(['r init =', num2str(sum_r_init)]);
    
    subplot(3,2,3)
    plot(0:15, avg_g_low.forward(1:16),'-og')
    hold on
    plot(0:-1:-15, avg_g_low.back(1:16),'-og');
    title(['g low =', num2str(sum_g_low)]);
    
    subplot(3,2,4)
    plot(0:15, avg_g_high.forward(1:16),'-og')
    hold on
    plot(0:-1:-15, avg_g_high.back(1:16),'-og');
    title(['g high =', num2str(sum_g_high)]);
    
    subplot(3,2,5)
    plot(0:15, avg_r_low.forward(1:16),'-or')
    hold on
    plot(0:-1:-15, avg_r_low.back(1:16),'-or');
    title(['r low =', num2str(sum_r_low)]);
    
    subplot(3,2,6)
    plot(0:15, avg_r_high.forward(1:16),'-or')
    hold on
    plot(0:-1:-15, avg_r_high.back(1:16),'-or');
    title(['r high =', num2str(sum_r_high)]);
    
    
    
    
    
    
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avg] = getAutoAvg(Auto)

mat_f = NaN*zeros(500,500);
mat_b = mat_f;

nn = numel(Auto);

    
for ii = 1:nn
    
    try 
    iif = numel(Auto{ii}.forward);
    iib = numel(Auto{ii}.back);
    
    mat_f(ii,1:iif) = Auto{ii}.forward;

    mat_b(ii,1:iib) = Auto{ii}.back;

    catch
    end
end
 

for ii = 1:500
    
    forward(ii) = mean(mat_f(~isnan(mat_f(:,ii)),ii));
    back(ii) = mean(mat_b(~isnan(mat_b(:,ii)),ii));

end

avg.forward = forward;
avg.back = back;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [sums] = corr_sums(tracks)

sums = sum(tracks.forward(~isnan(tracks.forward))) + ...
 sum(tracks.back(~isnan(tracks.back)));

end



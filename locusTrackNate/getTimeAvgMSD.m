function [MSD] = getTimeAvgMSD(dirname)

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

figure(20);
clf;
hold on;

figure(21);
clf;
hold on;

figure(22);
clf;
hold on;

for jj_list = 1:num_list_
    
    
    filename = [dirname,filesep,dir_list(jj_list).name];
    
    data = load(filename);
    
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
        
        msd.g_init = NaN;
        msd.r_init = NaN;
        msd.g_low = NaN;
        msd.g_high = NaN;
        msd.r_low = NaN;
        msd.r_high = NaN;
        
    else
        
        if isempty(sort.red.low) || isempty(sort.red.high)
            
            disp('Error: missing red track')
            
            msd.g_init = NaN;
            msd.r_init = NaN;
            msd.g_low = NaN;
            msd.g_high = NaN;
            msd.r_low = NaN;
            msd.r_high = NaN;
            
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
            
            %
            
            msd.g_init = getMSDint(g_init);
            msd.r_init = getMSDint(r_init);
            msd.g_low = getMSDint(g_low);
            msd.g_high = getMSDint(g_high);
            msd.r_low = getMSDint(r_low);
            msd.r_high = getMSDint(r_high);
            
            
            
            %         figure(20);
            %
            %         subplot(3,2,1)
            %         loglog(msd.g_init,'g');
            %         title('G initial')
            %         hold on
            %         subplot(3,2,2)
            %         loglog(msd.r_init,'r');
            %         title('R initial')
            %         hold on
            %         subplot(3,2,3)
            %         loglog(msd.g_low,'g');
            %         title('G Low')
            %         hold on
            %         subplot(3,2,4)
            %         loglog(msd.r_low,'r');
            %         title('R Low')
            %         hold on
            %         subplot(3,2,5)
            %         loglog(msd.g_high,'g');
            %         title('G High')
            %         hold on
            %         subplot(3,2,6)
            %         loglog(msd.r_high,'r');
            %         title('R High')
            %         hold on
            %
            
            x_int = log(msd.g_init);
            x_high = log(msd.g_high);
            t_ln = log(1:1:ss-2);
            
            if (ss <=100)
                F_h_ = polyfit(t_ln(~isnan(x_high(1:end)'))',x_high(~isnan(x_high(1:end)'))',1);
                
            else
                F_h_ = polyfit(t_ln(~isnan(x_high(1:80)'))',x_high(~isnan(x_high(1:80)'))',1);
            end
            
            F_i_ = polyfit(t_ln(~isnan(x_int(1:20)'))',x_int(~isnan(x_int(1:20)'))',1);
            
            F_h(jj_list) = F_h_(1);
            F_i(jj_list) = F_i_(1);
                       
            figure(21);
            %plot(t_ln,x_int,'g');
            loglog(msd.g_init,'-o','Color',[rand rand rand]);
            title('G initial');
            hold on;
            
            figure(22);
            %plot(t_ln,x_high,'g');
            loglog(msd.g_high,'-o','Color',[rand rand rand]);
            title('G High');
            hold on;

            
%             g_init(:,jj_list) = msd.g_init;
%             g_high(:,jj_list) = msd.g_high;
            
            
            
        end % END OF THE
    end     % SORT CHECK LOOPS
end

    MSD.high = F_h;
    MSD.int = F_i;
%    MSD.g_i = g_init;
%    MSD.g_h = g_high;

    figure(10);
    hist(F_h>0)
    
    figure(11);
    hist(F_i>0)

end

%%%%%%%%%%%%
%%%%%%%%%%%%

function [msd] = getMSDint(track)


ss = numel(track);


for itau = 1:(ss-2)
  
    for itime = 1:(ss-itau-2)
    
        r_(itime) = (track(itime+itau) - track(itime))^2;
    
    end
    
    ss_norm = 1/sum(~isnan(r_));
    msd(itau) = ss_norm*sum(r_(~isnan(r_)));
   
end


end



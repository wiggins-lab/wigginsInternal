function [window] = get_window(data_in)

%% THIS FUNCTION IS USED TO CALCULATE A RUNNING VARIANCE ON A TRAJECTORY IN
% AN EFFORT TO PICK OUT CHANGES IN DYNAMICS, I.E. A SHIFT FROM DIRECTED
% MOTION TO DIFFUSION

    sss = size(data_in);
    dat_var = [];

    
%% TRAVELLING LEFT TO RIGHT
    
    for nn = 1:sss(2)
    
        a = data_in(:,nn);
    
        a_ = a(~isnan(a));
    
        ss = numel(a_);

        box = 7;
    
        jj = 0;
    
        for iii =1:(ss-box)
        
            jjj = iii+box;
    
            dat_var_LR(iii,nn) = var(a_(iii:jjj));
    
        end

    
    end

    s = size(dat_var_LR);

    for ll = 1:s(1)
        
        dat_var_avg_LR(ll) = mean(dat_var_LR(ll,:));
    
    end
    
%% TRAVELLING RIGHT TO LEFT

    for nn = 1:sss(2)
    
        a = data_in(:,nn);
    
        a_ = a(~isnan(a));
    
        ss = numel(a_);

        box = 7;
    
        jj = 0;
    
        for iii =1:(ss-box-1)
        
            jjj = ss-iii;
            mmm = jjj-box;
    
            dat_var_RL(jjj,nn) = var(a_(mmm:jjj));
    
        end

    
    end

    s = size(dat_var_RL);

    for ll = 1:s(1)
        
        dat_var_avg_RL(ll) = mean(dat_var_RL(ll,:));
    
    end
    
%% %% METHOD 1 : Find first and last time variance goes above 5 %%%%    
    
%     t = 0:1:s(1);
%     
%     t = t(dat_var_avg >= 3); 
%     
%     y_minus = min(t);
%   
%     if (y_minus <= 0)
%         y_minus = 1;
%     end
%     
%     if ((max(t) + 3*box) <= s(1))
%         y_plus = max(t)+2*box ; 
%     else
%         y_plus = max(t) ;
%     end
%% %% METHOD 2 : Find max variance and expand around the peak    
%    
%     
%    [x,y] = max (dat_var_avg);
%         
%     y_minus = y ;
%    
%     if (y_minus < 0)
%         y_minus = 1;
%     end
%    
%     y_plus = y + 3*box;
%
%% %% METHOD 3 : Find max variance traveling in opposite directions

    [x_LR,y_LR] = max(dat_var_avg_LR);
    
    y_minus = y_LR;
    
    [x_RL,y_RL] = max(dat_var_avg_RL);

    y_plus = y_RL;

%% %% PLOT FOR DEBUGGING %%

    figure(20)
    clf
    plot(dat_var_avg_LR,'b')
    hold on
    plot(dat_var_avg_RL,'r')
    plot(y_minus,1:0.1:20,'-k')
    plot(y_plus,1:0.1:20,'-k')
    
   
    window.plus = y_plus;
    window.minus = y_minus;
    
end
function [slope] = getLocalSlope(track, window_)

window = window_+1;

ss = numel(track);

t = 1:ss;

ln = log10(track);
ln_t = log10(t);

%len = floor(ss/window);

 for ii = 1:(ss-window_)

    jj = ii + window_; 
     
    rise = ln(jj)-ln(ii);
    
    run = log10(window);
    
    slope(ii) = rise/run;
         
 end

slope = slope;

figure(10);
clf
plot(slope,'-o');

end





function [sort] = getTrackSort(xtrack)

% THIS FUNCTION IS DESIGNED TO SORT TRACKS IN XTRACK INTO UPPER AND LOWER
% TRACKS FOR CORRELATION ANALYSIS


low_ind = 0;
high_ind =0;

sort_low = [];
sort_high = [];

ss = size(xtrack);

for kk = 1:ss(2)
    
    a = xtrack((ss(1)-5):end,kk);
    avg = sum(a(~isnan(a)));
    
    if avg < 0       
        
        low_ind = low_ind + 1;
        sort_low(:,low_ind) = xtrack(:,kk);
        
    else
        
        high_ind = high_ind +1;
        sort_high(:,high_ind) = xtrack(:,kk);
        
    end
    
end



sort.low = sort_low;
sort.high = sort_high;

end
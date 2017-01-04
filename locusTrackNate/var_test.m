function [dat_varRL, dat_varLR] = var_test(track, box )


ss = numel(track);

box_ = floor(box/2);

for ii = (box_+1):(ss-box_)
    
    a = track((ii-box_):(ii+box_));
    
    a_ = a(~isnan(a));
    
    dat_varLR(ii) = var(a_);

%    t_var(ii) = t(ii);
end

%clf;
plot(dat_varLR,'r')


end
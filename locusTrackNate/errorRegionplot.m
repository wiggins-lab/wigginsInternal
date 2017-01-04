function errorRegionplot(msd,err,A,B)

%% THIS PROGRAM ADDS SHADED ERROR REGIONS TO A PLOT OF MEAN SQUARED
% DISPLACEMENT.  THE patch FUNCTION WILL NOT FILL IN THE REGION IF THERE
% ARE NaN'S IN THE ERROR LIST, SO YOU HAVE TO PASS THE REGION OF REAL
% NUMBERS (A TO B) WHICH THERE ARE NO NaN's IN THE err LIST
%
% NJK

t = 1:length(msd);

cc_ = [.8 .8 .8];
 
 patch([t(A:B) fliplr(t(A:B))], [msd(A:B) + ...
     err(A:B) fliplr(msd(A:B) - err(A:B))], cc_,'EdgeColor','none','FaceAlpha',0.5);


end
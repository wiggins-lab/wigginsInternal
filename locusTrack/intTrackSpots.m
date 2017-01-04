%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% intTrackSpots
%%
%% Fit locus positions from fluor image
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function locus = intTrackSpots(fluor, numLoci, celld, CONST);

opt =  optimset('MaxIter',25,'Display','off', 'TolX', 1/10);


% Does the dirty work of fitting loci
%locus = compSpotPosDev2(double(fluor),celld.mask,numLoci,crop,1,[],[],opt);
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,crop,1,[],[],opt);

%%
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,crop,1,[],[],opt);
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,CONST,1,[],[],opt);

%debug
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,CONST,1,1,2,opt);

locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,CONST,1,[],2,opt);

r=1;

%%

% correct the position of the loci with the offset of the lower corner of
% the image.
num_spot = length(locus);
for j = 1:num_spot;
    
    locus(j).r = locus(j).r + celld.r_offset-[1,1];
   
    locus(j).shortaxis = ...
        (locus(j).r-celld.coord.rcm)*celld.coord.e2;
    locus(j).longaxis = ...
        (locus(j).r-celld.coord.rcm)*celld.coord.e1;
    
end
end
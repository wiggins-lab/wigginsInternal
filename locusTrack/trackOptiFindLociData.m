%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% trackOptiMakeCell(dirname,numSpots)
%%
%% is part of the superSeggerOpti Package for tracking cells. This function
%% goes through the dirname/*err.mat files and computes the characteristics of the
%% cells in each frame and puts them in the CellA structure that is indexed
%% by the region number--not the cell ID. This code figures out the pole
%% age etc and also computes statistics on the fluorescence channels as
%% well as fitting the locus positions. These last two features are
%% controlled by parameters set in the loadConstants.m or
%% loadConstantsMine.m files.  numSpots variable in no longer used.
%%
%% Paul Wiggins, 4/17/2011
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_c = trackOptiFindLociData(data_c,CONST,header)

if ~exist('header')
    header = [];
end


nc = 0;
tmp_fn = fieldnames( data_c );
nf = numel( tmp_fn );
for j = 1:nf;
    if(strfind(tmp_fn{j},'fluor')==1)
        nc = nc+1;
    end
end

for ii = 1:data_c.regs.num_regs
    
    for j = 1:nc
        if isfield( CONST.trackLoci, 'numSpots' ) && numel(CONST.trackLoci.numSpots)>=j
            numSpots = CONST.trackLoci.numSpots(j);
            
            if numSpots
                tmp              = getfield( data_c.CellA{ii}, ['fluor',num2str(j)]);
                data_c.CellA{ii} = setfield( data_c.CellA{ii}, ['locus',num2str(j)], intTrackSpots( tmp, numSpots, data_c.CellA{ii}, CONST ));
            end
        end
    end
end

end
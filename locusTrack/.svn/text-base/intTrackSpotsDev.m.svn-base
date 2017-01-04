function celld = intTrackSpotsDev(celld,numLoci);

if (nargin< 2) || isempty(numLoci)
    NUM_SPOT_NUN = 2;
end

global CONST
if isempty( CONST )
    if exist('loadConstantsMine','file');
        loadConstantsMine
    else
        loadConstants
    end
end


crop = CONST.trackLoci.crop;


%
% Check to see if the fluor channel exist, then track loci.
%
if isfield(celld,'fluor1')
    if 0
        disp line
        celld.locus1 = compSpotPosDev2(double(celld.fluor1),celld.mask,numLoci,crop,1,1,7);
        drawnow;
        
     
        
        '';
    else
        % non-disp line
        celld.locus1 = compSpotPosDev2(double(celld.fluor1),celld.mask,numLoci,crop,1);
    end
end
% '';
if isfield(celld,'fluor2')
    if 0
        disp line
        celld.locus2 = compSpotPosDev2(double(celld.fluor2),celld.mask,numLoci,crop,1,1,7);
        drawnow;
        '';
    else
        % non-disp line
        celld.locus2 = compSpotPosDev2(double(celld.fluor2),celld.mask,numLoci,crop,1);
    end
end

if isfield(celld,'fluor3')
    celld.locus3 = compSpotPosDev(double(celld.fluo3),celld.mask,numLoci,crop,1);
end

if isfield(celld,'locus1');
    num_spot = length(celld.locus1);
    for j = 1:num_spot;
        celld.locus1(j).r = celld.locus1(j).r + celld.r_offset-[1,1];
        celld.locus1(j).shortaxis = ...
                (celld.locus1(j).r-celld.coord.rcm)*celld.coord.e2;
        celld.locus1(j).longaxis = ...
                (celld.locus1(j).r-celld.coord.rcm)*celld.coord.e1;
    end
end

if isfield(celld,'locus2');
    num_spot = length(celld.locus2);
    for j = 1:num_spot;
        celld.locus2(j).r = celld.locus2(j).r + celld.r_offset-[1,1];
        celld.locus2(j).shortaxis = ...
                (celld.locus2(j).r-celld.coord.rcm)*celld.coord.e2;
        celld.locus2(j).longaxis = ...
                (celld.locus2(j).r-celld.coord.rcm)*celld.coord.e1;    
    end
end

if isfield(celld,'locus3');
    num_spot = length(celld.locus3);
    for j = 1:num_spot;
        celld.locus3(j).r = celld.locus3(j).r + celld.r_offset(2:-1:1)-[1,1];
        celld.locus3(j).shortaxis=celld.coord.e1(1) * ( celld.locus3(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e1(2) * ( celld.locus3(j).r(1)-celld.coord.rcm(1) );
        celld.locus3(j).longaxis=celld.coord.e2(1) * ( celld.locus3(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e2(2) * ( celld.locus3(j).r(1)-celld.coord.rcm(1) );
    end
end


end
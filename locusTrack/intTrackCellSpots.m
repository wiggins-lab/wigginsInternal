function celld = intTrackCellSpots(celld);


if exist('loadConstantsMine','file');
    loadConstantsMine
else
    loadConstants
end



if R60X
    crop = 4;
elseif R100XB
    crop = 2;
else
    crop = 2;
end



if isfield(celld,'fluor1')
    if 0
        disp line
        celld.spot_1 = compSpotPos(double(celld.fluor1),celld.mask,crop,1,1,7);
        drawnow;
        '';
    else
        % non-disp line
        celld.spot_1 = compSpotPos(double(celld.fluor1),celld.mask,crop,1);
    end
end
% '';
if isfield(celld,'fluor2')
    if 0
        disp line
        celld.spot_2 = compSpotPos(double(celld.fluor2),celld.mask,crop,1,1,7);
        drawnow;
        '';
    else
        % non-disp line
        celld.spot_2 = compSpotPos(double(celld.fluor2),celld.mask,crop,1);
    end
end

if isfield(celld,'fluor3')
    celld.spot_3 = compSpotPos(double(celld.fluo3),celld.mask,crop,1);
end

if isfield(celld,'spot_1');
    num_spot = length(celld.spot_1);
    for j = 1:num_spot;
        celld.spot_1(j).r = celld.spot_1(j).r + celld.r_offset(2:-1:1)-[1,1];
        celld.spot_1(j).shortaxis=celld.coord.e1(1) * ( celld.spot_1(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e1(2) * ( celld.spot_1(j).r(1)-celld.coord.rcm(1) );
        celld.spot_1(j).longaxis=celld.coord.e2(1) * ( celld.spot_1(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e2(2) * ( celld.spot_1(j).r(1)-celld.coord.rcm(1) );
    end
end

if isfield(celld,'spot_2');
    num_spot = length(celld.spot_2);
    for j = 1:num_spot;
        celld.spot_2(j).r = celld.spot_2(j).r + celld.r_offset(2:-1:1)-[1,1];
        celld.spot_2(j).shortaxis=celld.coord.e1(1) * ( celld.spot_2(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e1(2) * ( celld.spot_2(j).r(1)-celld.coord.rcm(1) );
        celld.spot_2(j).longaxis=celld.coord.e2(1) * ( celld.spot_2(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e2(2) * ( celld.spot_2(j).r(1)-celld.coord.rcm(1) );
    end
end

if isfield(celld,'spot_3');
    num_spot = length(celld.spot_3);
    for j = 1:num_spot;
        celld.spot_3(j).r = celld.spot_3(j).r + celld.r_offset(2:-1:1)-[1,1];
        celld.spot_3(j).shortaxis=celld.coord.e1(1) * ( celld.spot_3(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e1(2) * ( celld.spot_3(j).r(1)-celld.coord.rcm(1) );
        celld.spot_3(j).longaxis=celld.coord.e2(1) * ( celld.spot_3(j).r(2)-celld.coord.rcm(2) )...
            +celld.coord.e2(2) * ( celld.spot_3(j).r(1)-celld.coord.rcm(1) );
    end
end


end
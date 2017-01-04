function [msd_L] = getCellLengthAvg( dirname )

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


for jj_list = 1:num_list_
    
    
    filename = [dirname,filesep,dir_list(jj_list).name];
    
    data = load(filename);
    
    [L,L0] = getCellLength( data );

    L_ = L/2 - L0/2;
    
    L_tot(1:numel(L_),jj_list) = L_.^2;
    
end

ss = size(L_tot);

for ii = 1:ss(1)
    
    msd_ = L_tot(ii,:);
    msd__ = msd_(msd_ > 0 );
    
    msd_L(ii) = mean(msd__);
    
end

end


function [L,L0] = getCellLength( data )

global CONST
PixelSize    = CONST.getLocusTracks.PixelSize;

ss = numel( data.CellA );

L0 = PixelSize*data.CellA{1}.length(1)*ones(1,ss);

for ii= 1:ss
    
    L(ii) = PixelSize*data.CellA{ii}.length(1);
    
end

end
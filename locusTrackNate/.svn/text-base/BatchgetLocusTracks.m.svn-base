function BatchgetLocusTracks( dirname , CONST)

%% THIS FUNCTION IS DESIGNED TO TAKE A DIRECTORY OF CELL FILES, FIND ALL
% THE TRACKS USING getLocusTracks_Dev.m, SAVE THE TRACKS TO THE DATA
% STRUCTURE, AND MOVE ALL THE FILES TO A NEW DIRECTORY CALLED tracks
%
% NOTE: CONVENIENT TO CALL THIS FROM A DIRECTORY ABOVE INSTEAD OF '.'
%
% NJK

if nargin < 2 || isempty('CONST')
    
    disp('Loading 60X Constants');
    CONST = loadConstants(60);

end

if ~exist('tracks','dir')
    
    mkdir tracks
    
end


dir_list = dir([dirname,filesep,'Cell*.mat']);

num_list_ = numel( dir_list );


for jj = 1:num_list_
    
    filename = [dirname,filesep,dir_list(jj).name]
    
    data = load( filename );
    
    if ~isfield(data,'track1') || ~isfield(data,'track2')
        
        disp([filename, ' : Getting New Locus Tracks']);
        dat = getLocusTracksDev_lite( data , [], CONST);
        
        save(['tracks', filesep, dir_list(jj).name], '-STRUCT' ,'dat');
        
    else
        
        disp([filename, ' : Using Existing Locus Tracks']);
        copyfile(filename,'tracks');
        
    end
    
end


end
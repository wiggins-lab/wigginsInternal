function fliptracks( filename )

% THIS FUNCTION FLIPS THE ORIENTATION OF THE TRACKS TO 
% CORRECT FOR 'LAZY' AND 'TRAVELLING' LOCI

data = load( filename );
data.track1.xtrack = -data.track1.xtrack;

if isfield(data,'locus2')
data.track2.xtrack = -data.track2.xtrack;
end

save( filename , '-STRUCT', 'data' );

end

% NOTE: THIS IS POSSIBLY A CODE TO GET Y POSITIONS
% test(:,2) = track_(ltrack(:,2));
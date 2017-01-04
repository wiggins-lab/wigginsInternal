
function [] = cutLoci( cell )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NJK, APRIL 2011 
% THIS FUNCTION IS AN INTERACTIVE TOOL TO EDIT LOCUS TRAJECTORIES 
% THAT ARE OUTPUT FROM 'multiLink'.  THIS PROGRAM ALLOWS YOU TO
% REMOVE ERRANT POINTS IN EACH TRAJECTORY.

% THINGS TO NOTE:
% - THIS FUNCTION WORKS ON Cell*.mat FILES THAT HAVE ALREADY
%       BEEN RUN THROUGH BatchgetLocusTracks or getLocusTracksDev
%
% - INPUT Cell*.mat file name  INTO THIS FUNCTION
%
% - EDITED TRACK FILES ARE SAVED IN THE SAME DIRECTORY, BUT WILL BE
%       OVERWRITTEN IF getLocusTracks IS RUN AGAIN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global CONST

if isempty( CONST )
    if exist('loadConstantsMine','file');
        loadConstantsMine
    else
        loadConstants
    end
end

%pixelsize = CONST.getLocusTracks.PixelSize;

pixelsize = 6/60;

filename = cell;

data = load(filename);

% CHECK TO SEE IF TRACKS HAVE ALREADY BEEN MADE
if ~isfield(data,'track1') || ~isfield(data,'track2')
    
    disp('NO TRACK DATA. TRY RUNNING BatchgetLocusTracks');
    
end

% CHECK DIRECTION OF POLE FOR PLOTTING PURPOSES
if exist('data.CellA{1}.pole.op_ori' )
    sign_ = sign( data.CellA{1}.pole.op_ori );
else
    sign_ = 1;
end

% DISPLAY THE CURRENT TRACKS OVER A KYMOGRAPH
kymoDisplay(data)

% CHOOSE WHICH TRACK TO EDIT
disp('Which track to edit?');
disp('g : Green');
disp('r : Red');

c = input(':', 's')

if c(1) == 'g'
    
    track = data.track1.xtrack;
    
else
    
    track = data.track2.xtrack;
    
end

% PICK ALL POINTS THAT YOU WANT TO REMOVE FROM SINGLE TRACK

disp('Choose consecutive points to remove');

x = ginput;

sss = size(x);

% DETERMINE THE TIME AND POSITIONS OF EACH POINT

for N = 1:sss(1);
    
    tind(N) = round(x(N,1));
    xind(N) = 10*x(N,2);
    
end

% FIND NUMBER OF TRACKS IN xtrack
% NOTE THAT MANY TRACKS HAVE SIMILAR POINTS THAT WE DON'T 
% WANT TO REMOVE FROM EVERY TRAJECTORY

ss = size(track);

% USE THE FIRST POINT TO DETERMINE WHICH TRACK TO EDIT

jj = 0;
% FOR EACH TRACK IN xtrack
for ii = 1:ss(2)
    
    % IF THE POSITION AT TIME N IS CLOSE TO THE OUTPUT OF ginput
    if abs(track(tind(1),ii) - xind(1)) < 1
        
        % SAVE THE COLUMN NUMBER (TRACK) THAT WE ARE CHANGING
        jj = jj+1;
        itrack(jj) = ii;
        
    end
    
end

% FOR EACH OF THE POINTS CHOSEN BEFORE
for nn = 1:sss(1)
    
    % FOR EACH OF THE TRACKS DETERMINED BY THE PREVIOUS LOOP
    for mm = 1:numel(itrack)
        
        % REPLACE ALL THE POSITIONS ASSOCIATED WITH THE TIMES
        % CHOSEN WITH NaNs
        track(tind(nn),itrack(mm)) = NaN;
        
    end
    
end

% PLOT EDITED TRACK ON KYMOGRAPH IN CYAN FOR CHECKING
figure(4)
xt = track;
ss = size(xt);
nt = ss(1);
ntracks = ss(2);
tt = 1:nt;
cc = 'c.';

for i = 1:ntracks
    xx = xt(:,i);    
    plot( sign_*pixelsize*xx, [cc,'-']);
    plot( tt(~isnan(xx)), sign_*pixelsize*xx(~isnan(xx)), [cc,':']);
end

% DECIDE WHETHER OR NOT TO KEEP THE NEW TRACK
disp('Keep new track data: Y or N ?')

d = input(':', 's')

if d(1) == 'y' || d(1) =='Y'
    
    % IF TRACK IS COPACETIC, REPLACE THE TRACK DATA IN XTRACK
    if c(1) == 'g'
        
        data.track1.xtrack = track;
        
    else
        
        data.track2.xtrack = track;
        
    end
    
    save( filename, '-STRUCT','data');
    
else
    
    % IF NO GOOD, DO NOTHING
    disp(['No Changes Made in ', filename]);

end

end


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

function kymoDisplay( data )
global pixelsize

num_frame = numel( data.CellA );

if exist('data.CellA{1}.pole.op_ori' )
    sign_ = sign( data.CellA{1}.pole.op_ori );
else
    sign_ = 1;
end

del = 1;

if ~isfield( data, 'kymo' );
    
    [Kymo,ll1] = makeKymograph( data, 0 );
    
    data.kymo = Kymo;
    data.kymo_ll = ll1;
    
else
    Kymo = data.kymo;
    ll1 = data.kymo_ll;
end

backer = autogain(Kymo.b);
backer = 0.15*(max(backer(:))-backer);

im = cat(3, autogain(Kymo.r )*del+backer, ...
    autogain(Kymo.g)*del+backer,...
    backer );

figure(4);
clf;
imagesc(1:num_frame, pixelsize*ll1, im );

track1 = data.track1;
track2 = data.track2;

figure(4)
hold on;
intPlotTrack( track1, 'g.' , sign_ )
intPlotTrack( track2, 'r.' , sign_ )

end

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

function intPlotTrack( track, cc, sign_ )

global pixelsize

xt = track.xtrack;
ss = size(xt);
nt = ss(1);
ntracks = ss(2);
tt = 1:nt;


for i = 1:ntracks
    xx = xt(:,i);
    
    plot( sign_*pixelsize*xx, [cc,'-']);
    plot( tt(~isnan(xx)), sign_*pixelsize*xx(~isnan(xx)), [cc,':']);
end


end

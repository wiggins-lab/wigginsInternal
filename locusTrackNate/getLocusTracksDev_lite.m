function [data] = getLocusTracks( data, kymo_disp, CONST );

%% THIS IS OUR BEST FUNCTION TO GET LOCI TRAJECTORIES.  THIS PARTICULAR
% VERSION IS SILIENT, I.E. DOES NOT OUTPUT ANY FIGURES DURING TRACK
% FINDING

%% INITIALIZE CONSTANTS

if isempty(CONST)
   
    disp('LOADING 60X CONSTANTS IN getLocusTracks');
    CONST = loadConstants(60);
    
end

FLUOR1_MIN_SCORE = CONST.getLocusTracks.FLUOR1_MIN_SCORE;
FLUOR1_REL       = CONST.getLocusTracks.FLUOR1_REL;

FLUOR2_MIN_SCORE = CONST.getLocusTracks.FLUOR2_MIN_SCORE;
FLUOR2_REL       = CONST.getLocusTracks.FLUOR1_REL;

pixelsize        = CONST.getLocusTracks.pixelsize;

if nargin < 2
    kymo_disp = 1;
end

num_frame = numel( data.CellA );

del = 1;

sign_ = sign( data.CellA{1}.pole.op_ori );

%% CHECK TO SEE IF TRACK DATA ALREADY EXISTS

if ~isfield(data,'track')
    
    n1 = 0;
    n2 = 0;
    %% get the size of the arrays;
    for i = 1:num_frame
        n1 = max( n1, numel(data.CellA{i}.locus1));
        
        if isfield(data,'locus2')
        n2 = max( n2, numel(data.CellA{i}.locus2));
        end
    end
    
    
    
    %% malloc arrays
    xx1 = zeros( num_frame, n1 )*NaN;
    yy1 = zeros( num_frame, n1 )*NaN;
    sc1 = zeros( num_frame, n1 )*NaN;
    fl1 = logical(zeros( num_frame, n1 ));
    
    
    if isfield(data,'locus2')
    xx2 = zeros( num_frame, n2 )*NaN;
    yy2 = zeros( num_frame, n2 )*NaN;
    sc2 = zeros( num_frame, n2 )*NaN;
    fl2 = logical(zeros( num_frame, n2 ));
    end
    
    ly = zeros( num_frame, 1 );
    lx = zeros( num_frame, 1 );
    
    
    %% Remove the tracks from the CellA structure
    for i = 1:num_frame
        for j = 1:numel(data.CellA{i}.locus1)
            
            xx1(i,j) = data.CellA{i}.locus1(j).longaxis;
            yy1(i,j) = data.CellA{i}.locus1(j).shortaxis;
            sc1(i,j) = data.CellA{i}.locus1(j).score;
        end
        fl1(i,:) = logical(double((sc1(i,:) > FLUOR1_MIN_SCORE))...
            + double(sc1(i,:) > FLUOR1_REL*max(sc1(i,:))));
        
        if isfield(data,'locus2')
        for j = 1:numel(data.CellA{i}.locus2)
            
            xx2(i,j) = data.CellA{i}.locus2(j).longaxis;
            yy2(i,j) = data.CellA{i}.locus2(j).shortaxis;
            sc2(i,j) = data.CellA{i}.locus2(j).score;
        end
        fl2(i,:) = logical(double((sc2(i,:) > FLUOR2_MIN_SCORE))...
            + double(sc2(i,:) > FLUOR2_REL*max(sc2(i,:))));
        end
        
        ly(i) = data.CellA{i}.length(2);
        lx(i) = data.CellA{i}.length(1);
    end
    
    data.track1.xx = xx1;
    data.track1.yy = yy1;
    data.track1.fl = fl1;
    data.track1.sc = sc1;
    data.track1.lx = lx;
    data.track1.ly = ly;
    
    if isfield(data,'locus2')
    data.track2.xx = xx2;
    data.track2.yy = yy2;
    data.track2.fl = fl2;
    data.track2.sc = sc2;
    data.track2.lx = lx;
    data.track2.ly = ly;
    end
else
    
    xx1  = data.track1.xx;
    yy1  = data.track1.yy;
    fl1  = data.track1.fl;
    sc1  = data.track1.sc;
    
    if isfield(data,'locus2')
    xx2  = data.track2.xx;
    yy2  = data.track2.yy;
    fl2  = data.track2.fl;
    sc2  = data.track2.sc;
    end
    
    lx   = data.track1.lx;
    ly   = data.track1.ly;
end

%% USE FUNCTION multiLink.m TO FIND TRAJECTORIES
track1       = data.track1;
track1.xx(~track1.fl) = NaN;
track1       = multiLink( track1 );
data.track1  = track1;
 
if isfield(data,'locus2')
track2       = data.track2;
track2.xx(~track2.fl) = NaN;
track2       = multiLink( track2 );
data.track2  = track2;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% intPlotTrack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% intPlotScore
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intPlotScore( track, cc )
global pixelsize

sc = track.sc;
fl = track.fl;

ss = size(sc);
nt = ss(1);
tt = 1:nt;
ntracks = ss(2);
tt = 1:nt;


for i = 1:ntracks
    scsc = sc(:,i);
    flfl = fl(:,i);
    
    plot( tt( flfl), scsc( flfl), [cc,'.-']);
    plot( tt(~flfl), scsc(~flfl), [cc,'o']);
    
end


end
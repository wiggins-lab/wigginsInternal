function [data] = getLocusTracks( data, kymo_disp, CONST );

if ~exist('CONST')
    disp('Loading Constants from loadConstants, 60X');
    CONST = loadConstants(60);
end


FLUOR1_MIN_SCORE = CONST.getLocusTracks.FLUOR1_MIN_SCORE;
FLUOR1_REL       = CONST.getLocusTracks.FLUOR1_REL;

FLUOR2_MIN_SCORE = CONST.getLocusTracks.FLUOR2_MIN_SCORE;
FLUOR2_REL       = CONST.getLocusTracks.FLUOR1_REL;

%pixelsize        = CONST.getLocusTracks.PixelSize;

pixelsize = 6/60;

if nargin < 2
    kymo_disp = 1;
end

%clf;

num_frame = numel( data.CellA );

del = 1;

if kymo_disp
    
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
    
    im1 = cat(3, 0*autogain(Kymo.r )*del+backer, ...
        autogain(Kymo.g)*del+backer,...
        backer );
    
    im2 = cat(3, autogain(Kymo.r )*del+backer, ...
        0*autogain(Kymo.g)*del+backer,...
        backer );
    
    figure(4);
    clf;
    imagesc(1:num_frame, pixelsize*ll1, im );
    
    figure(2);
    clf;
    imagesc(1:num_frame, pixelsize*ll1, im1 );
    
    
    figure(3);
    clf;
    imagesc(1:num_frame, pixelsize*ll1, im2 );
    
    
    figure(4);
    clf;
    imagesc(1:num_frame, pixelsize*ll1, im1 );
    
    figure(5);
    clf;
    imagesc((1:num_frame)*10/60, pixelsize*ll1, im2 );
    
    figure(6);
    clf;
    
    figure(7);
    clf
    
    
end


hold on;

if exist('data.CellA{1}.pole.op_ori' )
    sign_ = sign( data.CellA{1}.pole.op_ori );
else
    sign_ = 1;
end

%sign_ = 1;

if ~isfield(data,'track')
    
    n1 = 0;
    n2 = 0;
    % get the size of the arrays;
    for i = 1:num_frame
        n1 = max( n1, numel(data.CellA{i}.locus1));
        if isfield(data,'locus2')
            n2 = max( n2, numel(data.CellA{i}.locus2));
        end
    end
    
    
    
    % malloc arrays
    xx1 = zeros( num_frame, n1 )*NaN;
    yy1 = zeros( num_frame, n1 )*NaN;
    sc1 = zeros( num_frame, n1 )*NaN;
    fl1 = logical(zeros( num_frame, n1 ));
    
    xx2 = zeros( num_frame, n2 )*NaN;
    yy2 = zeros( num_frame, n2 )*NaN;
    sc2 = zeros( num_frame, n2 )*NaN;
    fl2 = logical(zeros( num_frame, n2 ));
    
    ly = zeros( num_frame, 1 );
    lx = zeros( num_frame, 1 );
    
    
    % Remove the tracks from the CellA structure
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
        end
        
        fl2(i,:) = logical(double((sc2(i,:) > FLUOR2_MIN_SCORE))...
            + double(sc2(i,:) > FLUOR2_REL*max(sc2(i,:))));
        
        
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
    
    xx2  = data.track2.xx;
    yy2  = data.track2.yy;
    fl2  = data.track2.fl;
    sc2  = data.track2.sc;
    
    lx   = data.track1.lx;
    ly   = data.track1.ly;
end


if kymo_disp
    %     figure(1)
    %     hold on;
    %
    %     plot(  lx/2*pixelsize, 'w')
    %     plot( -lx/2*pixelsize, 'w')
    %     plot( -lx/4*pixelsize, ':w')
    %     plot(  lx/4*pixelsize, ':w')
    %     plot(  lx*0*pixelsize, ':w')
    %
    %     xlabel( 'Frame (10 s)' )
    %     ylabel( 'Position (um)');
    
    %    figure(2)
    %    hold on;
    
    %     plot(  lx/2*pixelsize, 'w')
    %     plot( -lx/2*pixelsize, 'w')
    %     plot( -lx/4*pixelsize, ':w')
    %     plot(  lx/4*pixelsize, ':w')
    %     plot(  lx*0*pixelsize, ':w')
    %
    %     xlabel( 'Frame (10 s)' )
    %     ylabel( 'Position (um)');
    %
    %     figure(3)
    %     hold on;
    %
    %     plot(  lx/2*pixelsize, 'w')
    %     plot( -lx/2*pixelsize, 'w')
    %     plot( -lx/4*pixelsize, ':w')
    %     plot(  lx/4*pixelsize, ':w')
    %     plot(  lx*0*pixelsize, ':w')
    %
    %     xlabel( 'Frame (10 s)' )
    %     ylabel( 'Position (um)');
    
    
end





%figure(2)
%hold on;
disp('Finding Track 1:')
track1       = data.track1;
track1.xx(~track1.fl) = NaN;
track1       = multiLink( track1 );
%track1       = multiLinkRNA( track1 );

figure(4)
hold on;
track1.sumS = intPlotTrack( track1, [.4,1,.4] , sign_ );
data.track1  = track1;

if isfield(data,'locus2')
    hold on;
end
drawnow;


if isfield(data,'locus2')
    %figure(3)
    %hold on;
    disp('Finding Track 2:')
    track2       = data.track2;
    track2.xx(~track2.fl) = NaN;
    track2       = multiLink( track2 );
    track2.sumS = intPlotTrack( track2, [1,.4,.4] , sign_ );
    
    
    
    data.track2  = track2;
end



% figure(6);
% hold on;
% intPlotScore( track1, 'g' );
% intPlotScore( track2, 'r' );



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% intPlotTrack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sumS_ = intPlotTrack( track, cc, sign_ )
global pixelsize

xt = track.xtrack;
ss = size(xt);
nt = ss(1);
ntracks = ss(2);
tt = 1:nt;

meanS_ = [];
sumS_  = [];

for i = 1:ntracks
    xx = xt(:,i);
    
    plot( sign_*pixelsize*xx, ['.-'], 'Color', cc);
    plot( tt(~isnan(xx)), sign_*pixelsize*xx(~isnan(xx)), ['.:'],'Color', cc);
    
    
    if isfield( track, 'strack' )
        indEnd = find(  ~isnan(xx), 1, 'last' );
        
        sc = track.strack(:,i);
        
        meanS = mean(sc(~isnan(sc)));
        sumS = sum(sc(~isnan(sc)));
        text( indEnd+2, sign_*pixelsize*xx(indEnd), [num2str(meanS,'%2.2g'),', ',num2str(sumS,'%2.2g')], 'Color', 'w' );
        
        meanS_(i) = meanS;
        sumS_(i)  = sumS;
        
    end
    
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
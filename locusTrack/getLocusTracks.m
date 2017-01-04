function [data] = getLocusTracks( data, kymo_disp );

% 11-11-2010 NJK CHANGED `spot_1' and `spot_2' to 'locus1' and 'locus2'
% NOTE: There is another version of this program with 'spot' instead of locus
% called getLocusTracks_spot.m

global CONST

if exist('loadConstantsMine','file');
    loadConstantsMine(100)
    disp('LOADING')
else
    CONST = loadConstants( '100X', false );
end

TimeStep     = CONST.getLocusTracks.TimeStep;
pixelsize    = CONST.getLocusTracks.PixelSize;

%pixelsize = 6/60;

if nargin < 2
    kymo_disp = 1;
end

clf;

num_frame = numel( data.CellA );
time_step = TimeStep;
num_frame_step = num_frame*time_step;


%
% del sets the level of the background grey color.
%
del = .5;

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
    
    
    figure(1);
    clf;
    imagesc(1:time_step:num_frame_step, pixelsize*ll1, im );

    figure(4);
    clf;
    imagesc(1:time_step:num_frame_step, pixelsize*ll1, im1 );
    
    if exist('fluor2','var')
    
    figure(3);
    clf;
    imagesc(1:time_step:num_frame_step, pixelsize*ll1, im2 );

    end
end


hold on;


x11 = zeros( 1, num_frame )*NaN;
y11 = zeros( 1, num_frame )*NaN;
sc11= zeros( 1, num_frame )*NaN;
fl11= zeros( 1, num_frame )*NaN;

x12 = zeros( 1, num_frame )*NaN;
y12 = zeros( 1, num_frame )*NaN;
sc12= zeros( 1, num_frame )*NaN;
fl12= zeros( 1, num_frame )*NaN;

x21 = zeros( 1, num_frame )*NaN;
y21 = zeros( 1, num_frame )*NaN;
sc21= zeros( 1, num_frame )*NaN;
fl21= zeros( 1, num_frame )*NaN;

x22 = zeros( 1, num_frame )*NaN;
y22 = zeros( 1, num_frame )*NaN;
sc22= zeros( 1, num_frame )*NaN;
fl22= zeros( 1, num_frame )*NaN;

ly = zeros( 1, num_frame );
lx = zeros( 1, num_frame );


% Copy the tracks from the CellA structure
% into arrays that can be plotted
for i = 1:num_frame
    if ~isempty( data.CellA{i}.locus1 )
        
        x11(i) = data.CellA{i}.locus1(1).longaxis;
        y11(i) = data.CellA{i}.locus1(1).shortaxis;
        sc11(i)= data.CellA{i}.locus1(1).score;
    else
        x11(i) = NaN;
        y11(i) = NaN;
        sc11(i)= NaN;
    end
    
    if numel( data.CellA{i}.locus1 ) > 1
        x12(i) = data.CellA{i}.locus1(2).longaxis;
        y12(i) = data.CellA{i}.locus1(2).shortaxis;
        sc12(i)= data.CellA{i}.locus1(2).score;
    else
        x12(i) = NaN;
        y12(i) = NaN;
        sc12(i)= NaN;
        
    end
    
    if exist('locus2','var')
    
    if ~isempty( data.CellA{i}.locus2 )
        x21(i) = data.CellA{i}.locus2(1).longaxis;
        y21(i) = data.CellA{i}.locus2(1).shortaxis;
        sc21(i)= data.CellA{i}.locus2(1).score;
    else
        x21(i) = NaN;
        y21(i) = NaN;
        sc21(i)= NaN;
    end
    
    
    if numel( data.CellA{i}.locus2 ) > 1
        x22(i) = data.CellA{i}.locus2(2).longaxis;
        y22(i) = data.CellA{i}.locus2(2).shortaxis;
        sc22(i)= data.CellA{i}.locus2(2).score;
    else
        x22(i) = NaN;
        y22(i) = NaN;
        sc22(i)= NaN;
    end
    
    end
    
    ly(i) = data.CellA{i}.length(2);
    lx(i) = data.CellA{i}.length(1);
end



track.lx = lx;
track.ly = ly;

% the x var is the long axis position
xx1 = [x11',x12'];
xx2 = [x21',x22'];

% the y var is the short axis position
yy1 = [y11',y12'];
yy2 = [y21',y22'];

% sc is the score.
sc1 = [sc11',sc12'];
sc2 = [sc21',sc22'];


% make sure that the new pole is in the + direction.
%

sign_ = sign( data.CellA{1}.pole.op_ori );
if sign_ == 0
    sign_ = 1;
end


% set plot flags based on limits set in
% 11-11-2010 CONSTANTS FOR getLocusTracks DON'T SEEM TOO LOAD 
% THROUGH loadConstantsMine SO ADDED THEM EXPLICITLY HERE
FLUOR1_MIN_SCORE = 2;
FLUOR2_MIN_SCORE = 2;
FLUOR1_REL = 0.3;
FLUOR2_REL = 0.3;

fl11 = sc11 > FLUOR1_MIN_SCORE;
fl12 = logical(double((sc12 > FLUOR1_MIN_SCORE)) .* double(sc12 > FLUOR1_REL*sc11));

fl21 = sc21 > FLUOR2_MIN_SCORE;
fl22 = logical(double((sc22 > FLUOR2_MIN_SCORE)) .* double(sc22 > FLUOR2_REL*sc21));

fl1 = [fl11',fl12'];
fl2 = [fl21',fl22'];

% need to fix switching
% This function peices together the individual points into locus tracks
[map1] = fixTrack( xx1 );
xx1 = mapIt(xx1, map1);
yy1 = mapIt(yy1, map1);
sc1 = mapIt(sc1, map1);
fl1 = mapIt(fl1, map1);

[map2] = fixTrack( xx2 );
xx2 = mapIt(xx2, map2);
yy2 = mapIt(yy2, map2);
sc2 = mapIt(sc2, map2);
fl2 = mapIt(fl2, map2);

% Plot trcks and Kymo.

t = 1:numel(xx1(:,1));
t_step = time_step*t;
figure(1);
hold on;
if kymo_disp
     plot( t_step, lx/2*pixelsize, ':w')
     hold on;
     plot( t_step, -lx/2*pixelsize, ':w')
     plot( t_step, -lx/4*pixelsize, ':w')
     plot( t_step, lx/4*pixelsize, ':w')
     plot( t_step, lx*0*pixelsize, ':w')
    
     plot( t_step, sign_*xx1(:,1)*pixelsize, 'g:')
     plot( t_step,sign_*xx2(:,1)*pixelsize, 'r:')
     plot( t_step,sign_*xx1(:,2)*pixelsize, 'c:')
     plot( t_step,sign_*xx2(:,2)*pixelsize, 'm:')
     
    xx1_ = xx1;
    xx1_(~fl1(:,1),1) = NaN;
    xx1_(~fl1(:,2),2) = NaN;
    
    xx2_ = xx2;
    xx2_(~fl2(:,1),1) = NaN;
    xx2_(~fl2(:,2),2) = NaN;
%     
    
figure(4)
hold on;
     plot( t_step,sign_*xx1_(:,1)*pixelsize, 'g.')
%     plot( t,xx2_(:,1)*pixelsize, 'r-')
     plot( t_step,sign_*xx1_(:,2)*pixelsize, 'g.')
%     plot( t,xx2_(:,2)*pixelsize, 'm-')
    
     x2x2 = -abs(xx1_(:,2) - xx1_(:,1));
   % plot( t_step, x2x2*pixelsize, '-w.');
    
    xlabel( 'Time (min)' )
    ylabel( 'Position (um)');

    %end


% ADDED BY NJK 11-16-2010 to put red dots on M-Cherry graph
% and RED and GREEN dots on combined graph
% 11-30-10 THE Y-AXIS (X POSITIONS) ARE REVERSED

if exist('locus2','var')
figure(3)
hold on;
     plot( t_step,sign_*xx2_(:,1)*pixelsize, 'r.')
%     plot( t,xx2_(:,1)*pixelsize, 'r-')
     plot( t_step,sign_*xx2_(:,2)*pixelsize, 'r.')
%     plot( t,xx2_(:,2)*pixelsize, 'm-')
    
    x2x2 = -abs(xx2_(:,2) - xx2_(:,1));
   % plot( t_step, x2x2*pixelsize, '-w.');

    
    xlabel( 'Time (min)' )
    ylabel( 'Position (um)');

end
    
figure(1)
hold on;
     plot( t_step,sign_*xx2_(:,1)*pixelsize, 'r.')
%     plot( t,xx2_(:,1)*pixelsize, 'r-')
     plot( t_step,sign_*xx2_(:,2)*pixelsize, 'r.')
%     plot( t,xx2_(:,2)*pixelsize, 'm-')
     plot( t_step,sign_*xx1_(:,1)*pixelsize, 'g.')
%     plot( t,xx2_(:,1)*pixelsize, 'r-')
     plot( t_step,sign_*xx1_(:,2)*pixelsize, 'g.')    
    
    xlabel( 'Time (min)' )
    ylabel( 'Position (um)');

end

% figure(2)
% clf;
% plot(sc1(:,1), 'g-')
% hold on;
% plot(sc1(:,2), 'c-')
% plot(sc2(:,1), 'r-')
% plot(sc2(:,2), 'm-')


track = [];

track.x1  = xx1;
track.y1  = yy1;
% track.sc1 = sc1;
% track.fl1 = fl1;
% 
% 
track.x2  = xx2;
track.y2  = yy2;
% track.sc2 = sc2;
% track.fl2 = fl2;

% x1 is the green track

track.x11  = sign_*xx1_(:,1)*pixelsize;
track.x12  = sign_*xx1_(:,2)*pixelsize;
track.y1  = yy1;
track.sc1 = sc1;
track.fl1 = fl1;

% x2 is the red track

track.x21  = sign_*xx2_(:,1)*pixelsize;
track.x22  = sign_*xx2_(:,2)*pixelsize;
track.y2  = yy2;
track.sc2 = sc2;
track.fl2 = fl2;



track.time = t_step;

data.track = track;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Connect tracks into continuous peices...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map] = fixTrack( xx );
map = xx*0;


dist_ = [(xx(2:end,:)-xx(1:end-1,:)).^2];
distX = [(xx(2:end,1)-xx(1:end-1,2)).^2,(xx(2:end,2)-xx(1:end-1,1)).^2];

dist= [sum(dist_')',sum(distX')'];

ns = size(xx);

nn = ns(1);
sflag = 0;

map(1,:) = [1,2];

for i = 2:nn
    if all(~isnan(dist_(i-1,:)))
        if dist(i-1,1) > dist(i-1,2)
            sflag = ~sflag;
        end
    elseif  isnan(xx(i-1,1)) && ~isnan(xx(i-1,2))...
            && ~isnan(xx(i,1)) && ~isnan(xx(i,2))
        if (xx(i-1,2)-xx(i,1))^2 < (xx(i-1,2)-xx(i,2))^2
            sflag = ~sflag;
        end
    elseif ~isnan(xx(i-1,1)) &&  isnan(xx(i-1,2))...
            && ~isnan(xx(i,1)) && ~isnan(xx(i,2))
        if (xx(i-1,1)-xx(i,1))^2 > (xx(i-1,1)-xx(i,2))^2
            sflag = ~sflag;
        end
    elseif ~isnan(xx(i-1,1)) && ~isnan(xx(i-1,2))...
            &&  isnan(xx(i,1)) && ~isnan(xx(i,2))
        if (xx(i,2)-xx(i-1,2))^2 > (xx(i,2)-xx(i-1,1))^2
            sflag = ~sflag;
        end
    elseif ~isnan(xx(i-1,1)) && ~isnan(xx(i-1,2))...
            && ~isnan(xx(i,1)) &&  isnan(xx(i,2))
        if (xx(i,1)-xx(i-1,1))^2 > (xx(i,1)-xx(i-1,2))^2
            sflag = ~sflag;
        end
    else
        %'huh'
    end
    
    if sflag
        map(i,:) = [2,1];
    else
        map(i,:) = [1,2];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mapIt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx = mapIt( xx, map )

ss = size(xx);

for i = 1:ss(1)
    xx(i,:) = xx(i,map(i,:));
end
end
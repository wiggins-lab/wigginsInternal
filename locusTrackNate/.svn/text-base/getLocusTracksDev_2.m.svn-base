function [data] = getLocusTracksDev_2( dirname, CONST );

%% THIS PROGRAM IS USED TO CREATE A LARGE FIGURE WITH KYMOGRAPHS FOR ALL
% COMPLETE CELL CYCLES IN dirname PLOTTED WITH THE TRAJECTORIES OVERLAID.
% IT PRODUCES A BUNCH OF KYMOGAPHS, STORES THEM IN A TEMP DIRECTORY, THEN
% RE-READS THEM IN AND ULTIMATELY REMOVES THE TEMP DIRECTORY
% 
% NJK

%% INITIALIZE CONSTANTS

if isempty( CONST )
    
    disp('LOADING 60X CONSTANTS IN getLocusTracks');
    CONST = loadConstants(60);
    
end

%% Define parameters

FLUOR1_MIN_SCORE = CONST.getLocusTracks.FLUOR1_MIN_SCORE;
FLUOR1_REL       = CONST.getLocusTracks.FLUOR1_REL;

FLUOR2_MIN_SCORE = CONST.getLocusTracks.FLUOR2_MIN_SCORE;
FLUOR2_REL       = CONST.getLocusTracks.FLUOR1_REL;

pixelsize        = CONST.getLocusTracks.pixelsize;


%% Read in files

dir_list = dir([dirname,filesep,'Cell*.mat']);

mkdir tmp;

num_list_ = numel( dir_list );
num_list = 0;

for ii = 1:num_list_
    
    if (~isempty(dir_list(ii).name)) && ( dir_list(ii).name(1) == 'C' )
        num_list = num_list + 1;
    end
end

numa = ceil(sqrt(num_list));
numb = ceil( num_list/numa);

if nargin < 3
    kymo_disp = 1;
end

%% Begin Loop
for jj = 1:num_list_
    
    if (~isempty(dir_list(jj).name)) && ( dir_list(jj).name(1) == 'C' )
        ii = ii + 1;
        
        
        filename = [dirname,filesep,dir_list(jj).name]
        
        data = load(filename);
        
        %% GET TRACK STRUCTURE, INCLUDING XTRACK
        %   FIRST CHECK TO SEE IF TRACKS HAVE ALREADY BEEN CREATED
        
        if ~isfield(data,'track1') || ~isfield(data,'track2')
            
            % OPTION 1: MAKE NEW TRACKS FOR EACH KYMO
            
            disp([filename, ' : Getting New Locus Tracks']);
            
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
                
                figure(4);
                clf;
                imagesc(1:num_frame, pixelsize*ll1, im );
                
                
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
                
                
                
                data.track2.xx = xx2;
                data.track2.yy = yy2;
                data.track2.fl = fl2;
                data.track2.sc = sc2;
                data.track2.lx = lx;
                data.track2.ly = ly;
            else
                
                xx1  = data.track1.xx;
                yy1  = data.track1.yy;
                fl1  = data.track1.fl;
                sc1  = data.track1.sc;
                
                if isfield(data.CellA{1},'locus2')
                xx2  = data.track2.xx;
                yy2  = data.track2.yy;
                fl2  = data.track2.fl;
                sc2  = data.track2.sc;
                end
                
                lx   = data.track1.lx;
                ly   = data.track1.ly;
            end
            
            
            
            % figure(2)
            % hold on;
            % disp('Process 1:')
            track1       = data.track1;
            track1.xx(~track1.fl) = NaN;
            track1       = multiLink( track1 );
            data.track1  = track1;
            %
            if isfield(data.CellA{1},'locus2')
            % %if exist('locus2','var')
            % figure(3)
            % hold on;
            % disp('Process 2:')
            track2       = data.track2;
            track2.xx(~track2.fl) = NaN;
            track2       = multiLink( track2 );
            data.track2  = track2;
            end
            
            figure(4)
            hold on;
            intPlotTrack( track1, 'g.' , sign_ )
            if isfield(data.CellA{1},'locus2')
            intPlotTrack( track2, 'r.' , sign_ )
            end
            title(['Cell ',num2str(data.ID)]);

            %set(gcf,'Color',[0 0 0]);
            
        else
            
            
            % OPTION 2: USE EXISTING TRACKS
            
            disp([filename, ' : Using Existing Locus Tracks']);
            kymoDisplay(data);
            title(['Cell ',num2str(data.ID)]);
            
            %set(gcf,'Color',[0 0 0]);
        end
        
        
        
        
        print(4, '-dtiff', ['tmp/', num2str(jj), '.tiff'])
        
        cell_id(jj) = data.ID;
        
    end
    
end

%% PLOT ALL KYMOGRAPHS ON A SINGLE FIGURE

figure(10);
clf;

plot_num = ceil(sqrt(num_list_));
plot_num2 = plot_num*plot_num;




for kk = 1:plot_num2
    
    im_combo{kk} = zeros(900,1200, 3);
    
end

for kk = 1:num_list_
    
    im_combo{kk} = imread(['tmp/',num2str(kk),'.tiff']);
    
end
 try 
for ll = 1:plot_num
    
    im_{ll} = [];
    
end
 catch
     keyboard
 end

for mm = 1:plot_num2
    
    jj = floor((mm-1)/(plot_num)) + 1;
    
    im_{jj} = cat(2, im_{jj}, im_combo{mm});
    
end



im_tot = [];

for nn = 1:plot_num
    
    im_tot = cat(1,im_tot,im_{nn});
    
end

imshow(im_tot);
%text(0,0,num2str(cell_id(kk)));

%% REMOVE THE TEMPORY DIRECTORY WHERE IMAGE FILES WERE STORED

!rm -r tmp/


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
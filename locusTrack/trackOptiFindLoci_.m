%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% trackOptiMakeCell(dirname,numSpots)
%%
%% is part of the superSeggerOpti Package for tracking cells. This function
%% goes through the dirname/*err.mat files and computes the characteristics of the
%% cells in each frame and puts them in the CellA structure that is indexed
%% by the region number--not the cell ID. This code figures out the pole
%% age etc and also computes statistics on the fluorescence channels as
%% well as fitting the locus positions. These last two features are
%% controlled by parameters set in the loadConstants.m or
%% loadConstantsMine.m files.  numSpots variable in no longer used.
%%
%% Paul Wiggins, 4/17/2011
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [clist, clist_def] = trackOptiFindLoci(dirname,CONST,header)

if ~exist('header')
    header = [];
end

dirname = fixDir( dirname );

% Get the track file names...
contents=dir([dirname  '*_err.mat']);

num_im = numel(contents);

data_c = loaderInternal([dirname  contents(1  ).name]);

nc = 0;
tmp_fn = fieldnames( data_c );
nf = numel( tmp_fn );
for j = 1:nf;
    if(strfind(tmp_fn{j},'fluor')==1)
        nc = nc+1;
    end
end

% Find Loci can be ||
if isfield( data_c, 'fluor1' );
    if CONST.show_status
        h = waitbar( 0, 'Find Loci.');
    else
        h = [];
    end
    for i = 1:num_im
        %parfor i = 1:num_im;
        intDoLoci( i, dirname, contents, num_im, nc, CONST);
        if CONST.show_status
            waitbar(i/num_im,h,['Find Loci--Frame: ',num2str(i),'/',num2str(num_im)]);
        else
            disp( [header, 'FindLoci: No status bar. Frame ',num2str(i), ...
                ' of ', num2str(num_im),'.']);
        end
    end
    if CONST.show_status
        close(h);
    end
end

end

function data = loaderInternal( filename );
data = load( filename );
end




function intDoLoci( i, dirname, contents, num_im, nc, CONST)

data_c = loaderInternal([dirname,contents(i  ).name]);
gf  = fspecial( 'gaussian', 21, 3 );
gf2 = fspecial( 'gaussian', 21, 2 );

% make filtered images to fit
im_filt = cell([1,nc]);
for j = 1:nc
    fl_im = getfield( data_c, ['fluor',num2str(j)]);
    fl_im = medfilt2( double(fl_im), [3,3], 'symmetric' );
    
    Istd  = std(fl_im(data_c.mask_bg));
    tmp = (fl_im-imfilter( fl_im, gf, 'replicate' ))/Istd;
    tmp(tmp<0) = 0;
    im_filt{j} = tmp;
end

backer =  imfilter( sqrt(bwdist( ~data_c.mask_cell )), gf2, 'replicate' );


for ii = 1:data_c.regs.num_regs
    
    
%     if ii == 11
%    'hi'    
%     end
    
    for j = 1:nc
        if isfield( CONST.trackLoci, 'numSpots' ) && numel(CONST.trackLoci.numSpots)>=j
            numSpots = CONST.trackLoci.numSpots(j);
            
            if numSpots
                fluor = getfield( data_c.CellA{ii}, ['fluor',num2str(j)]);
                xx    = data_c.CellA{ii}.xx;            
                yy    = data_c.CellA{ii}.yy;            

                ff = im_filt{j}(yy,xx);
                backer_ = backer(yy,xx);
                
                data_c.CellA{ii} = setfield( data_c.CellA{ii}, ['locus',num2str(j)], ...
                    intTrackSpots2( fluor, ff, backer_, numSpots, data_c.CellA{ii}, CONST ));
            end
        end
    end
end

dataname = [dirname,contents(i  ).name];
save(dataname,'-STRUCT','data_c');

end

function locus = intTrackSpots2(fluor, ff, backer, numLoci, celld, CONST);

opt =  optimset('MaxIter',25,'Display','off', 'TolX', 1/10);


% Does the dirty work of fitting loci
%locus = compSpotPosDev2(double(fluor),celld.mask,numLoci,crop,1,[],[],opt);
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,crop,1,[],[],opt);

%%
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,crop,1,[],[],opt);
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,CONST,1,[],[],opt);

%debug
%locus = compSpotPosDevFmin(double(fluor),celld.mask,numLoci,CONST,1,1,2,opt);
disp_flag = false;

locus = compSpotPosDevFmin3(double(fluor), double(ff), double(backer), ...
    celld.mask, numLoci,CONST,1,disp_flag,2,opt);

r=1;

%%

% correct the position of the loci with the offset of the lower corner of
% the image.
num_spot = length(locus);
for j = 1:num_spot;
    
    locus(j).r = locus(j).r + celld.r_offset-[1,1];
   
    locus(j).shortaxis = ...
        (locus(j).r-celld.coord.rcm)*celld.coord.e2;
    locus(j).longaxis = ...
        (locus(j).r-celld.coord.rcm)*celld.coord.e1;
    
end
end


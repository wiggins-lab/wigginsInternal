%% Model of cytoplasmic fluor in cell. Fit cell by cell.
% numc is channel number
% function adds the field cytoX, where cyto is the global cytofluor model
% for channel X = numc.
function [ data ] = intFindFociCyto( data, CONST, numc )

disp_flag = false;
disp_flag2 = false;


% Get images out of the structures.
image0 = double(getfield( data, ['fluor',num2str(numc)] ));
cyto0 = double(getfield( data, ['cyto',num2str(numc)] ));
image = double(image0)-double(cyto0);
image(image<0) = 0;

% kill cyto stuff
%image = image0;
%image = medfilt2( image0, [3,3], 'symmetric' );


dI0 = std( image(data.mask_bg) );

% make a medfilt highpass version of the image.
gf_hp  = fspecial( 'gaussian', 21, 3 );
image0_med = medfilt2( image, [3,3], 'symmetric' );

Istd  = std(image(data.mask_bg));
image0_medhp = (image0_med-imfilter( image0_med, gf_hp, 'replicate' ))/Istd;
image_fil = image0_medhp;
image_fil(image_fil<0) = 0;

% show the image
% imshow( image_fil, [] );

if disp_flag2
    figure(1);
    clf;
    imshow( image0, [] );
    
    hold on;
    
end


%% Stroke out the mask by two pixel
mask_mod = bwmorph( data.mask_bg, 'dilate', 2 );


%% Set up the persistent variables here
fieldname = ['locus',num2str(numc)];

gaussR = CONST.trackLoci.gaussR;
crop   = CONST.trackLoci.crop;

opt =  optimset('MaxIter',25,'Display','off', 'TolX', 1/10);


focus0 = struct(...
    'r',               [NaN,NaN], ...
    'score',           NaN, ...
    'error',           NaN, ...
    'intensity_score', NaN, ...
    'intensity',       NaN, ...
    'b',               NaN,...
    'shortaxis',       NaN,...
    'longaxis',        NaN );






sdisk =  strel('disk',1);
gauss  =  fspecial('gaussian',7,gaussR);
%lap    = -fspecial('laplacian');
%locusd =  imfilter(gauss, lap);
crop_ = -crop:crop;
[Xpre, Ypre] = meshgrid(crop_,crop_);



%% Make padded image.
% pad the images to facilitate cropping around each spot
%
% there are two masks: one for watersheding and one for fitting.
% the fitting mask is slightly wider to accomodate spots near
% the cell boundary.

% input image size
imsize     = size(image);
% expanded image size
imsize_    = imsize +  crop*2;


if ~exist('disp_flag', 'var');
    disp_flag = 0;
end
if disp_flag
    if ~exist('fig_num','var')
        fig_num = 7;
    end
end


% expanded image
im_    = zeros( imsize_ );
im_(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = image;

model_ = im_;

im0_    = zeros( imsize_ );
im0_(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = image0;

imf_    = zeros( imsize_ );
imf_(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = image_fil;

% expanded mask
mk_     = false( imsize_ );
mk_(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = mask_mod;

% make a smoothed image to watershed
im0_s_ = imfilter(im0_, gauss, 'replicate');

% zero out regions of fluor outside of cells.
im0_s_(~mk_)=0;

%%%%%%%%
im0_s_ = im0_s_-dI0;
im0_s_(im0_s_<0) = 0;
%%%%%%%%%

% perform the wateshed
L = watershed(-im0_s_);
L(~mk_) = 0;
mask__ = logical(L);
%mask_ = bwmorph( bwmorph( bwmorph( mask__, 'fill' ), 'erode', 2 ), 'dilate', 2 );
mask_ =  bwmorph( mask__, 'fill' );

%% Striop out small shit
fluor_regs = bwlabel( mask_ );
fluor_props = regionprops( fluor_regs, {'Area', 'BoundingBox'} );


dx = drill(  fluor_props, '.BoundingBox(3)' );
dy = drill(  fluor_props, '.BoundingBox(4)' );
A  = drill(  fluor_props, '.Area' );


AreaLowerBound  = 10;
WidthLowerBound = 4;

keepers = find(all( [dx>=WidthLowerBound,dy>=WidthLowerBound,A>=AreaLowerBound],2 ));

mask_keeps = ismember(fluor_regs, keepers );

fluor_regs = bwlabel(mask_keeps );
fluor_props = regionprops( fluor_regs, imf_, {'Area', 'BoundingBox', 'MaxIntensity'} );
num_regs = max(fluor_regs(:));


% Loop through the regions and fit loci in each one
for ii = 1:num_regs
    
    
    [xx,yy] = getBBpad( fluor_props(ii).BoundingBox, imsize_, crop );
    
    
    im0_ii  = im0_(yy,xx);
    im_ii   = im_(yy,xx);
    imf_ii   = imf_(yy,xx);
    
    ss_ii = size( im0_ii );
    
    mask_ii = (fluor_regs(yy,xx)==ii);
    
    im_tmp = imf_ii;
    im_tmp( ~mask_ii ) = nan;
    
    [maxValue,y] = max(im_tmp(:));
    
    [y,x] = ind2sub(ss_ii,y);
    
    xn = x;
    yn = y;
    fn = maxValue;
    
    an = sum(mask_ii(:));
    
    
    %%
    %---------------------------------------------------------
    %
    % plots for debugging
    %
    %---------------------------------------------------------
    if disp_flag
        figure(5);
        clf;
        
        backer = ag( im0_ii );
        
        imshow( cat( 3, ...
            0.33* ag( ~mask_ii ) + backer, ...
            backer, ...
            backer) );
        
        hold on;
        
        % plot guess
        plot( x, y, 'r*' );
        
    end
    
    
    
    % now we fit a gaussian to locate the spot within each region
    
    % Initial values loaded in here
    x0 = xn;
    y0 = yn;
    I0 = 0;
    
    IG0 = im_ii(y,x)-I0;
    b0 = 1;
    %bmin = 1;
    %bmax = 6;
    
    ddd = [x0,y0,IG0,b0];
    
    % lower bound
    %lddd = [x-crop,y-crop,0,bmin];
    
    % upper bound
    %uddd = [x+crop,y+crop,100000,bmax];
    
    % make the coordinates for the cropped image
    ix2 = x+crop_;
    iy2 = y+crop_;
    
    xx2 = Xpre+x;
    yy2 = Ypre+y;
    
    
    % here we actually do the fit
    
    
    mkk = mask_ii(iy2,ix2);
    imm = im_ii(iy2,ix2);
    
    
    
    
    [ddd,res] = fminsearch( @do_the_fit,ddd,opt);
    
    
    % extract the position, width, and intensity of the fitted gaussian
    
    x  = ddd(1)+xx(1)-1-crop;
    y  = ddd(2)+yy(1)-1-crop;
    IG = ddd(3);
    b  = ddd(4);
    
    
    % subtract out the focus from the model
    [xx__,yy__] = meshgrid( xx, yy );
    model__ = g_fit_fun(xx__,yy__,x+crop,y+crop,ddd(3),I0,ddd(4));
    model_(yy,xx) = model_(yy,xx) - model__;
    
    
    if disp_flag
        figure(77);
        clf;
        imshow( [model__, im_ii, model_(yy,xx)], [] );
    end
    
    
    
    bs = .5;
    score = IG/dI0; %*exp( -(b-1).^2/(2*bs^2));
    
    focus = focus0;
    
    focus.r               = [x,y];
    
    focus.score           = score;
    focus.intensity_score = IG/dI0;
    focus.intensity       = IG;
    focus.b               = b;
    focus.error           = sqrt(res/sum(mkk(:))/IG.^2);
    
    
    %% figure out which cell
    x_  = min( [max( [1,round( x )]),imsize(2)] );
    xmin = max( [1,round( x )-5] );
    xmax = min( [imsize(2),round( x )+5]);
    if xmin>xmax;
        xmin = xmax;
    end
    
    y_  = min( [max( [1,round( y )]),imsize(1)] );
    ymin = max( [1,round( y )-5] );
    ymax = min( [imsize(1),round( y )+5]);
    if ymin>ymax;
        ymin = ymax;
    end
    
    xx = xmin:xmax;
    yy = ymin:ymax;
    
    ss_dist = [numel(yy),numel(xx)];
    
    cells_label = data.regs.regs_label(yy,xx);
    cells_mask = logical(cells_label);
    mask_tmp = zeros( ss_dist );
    mask_tmp( y_-yy(1)+1, x_-xx(1)+1 ) = 1;
    dist_tmp = bwdist( mask_tmp );
    
    list = cells_label(cells_mask);
    [~,ind_min] = min( dist_tmp(cells_mask) );
    cell_num = list(ind_min);
    
    %% Add the focus to the right cell
    if ~isempty( cell_num )
        if isfield( data.CellA{cell_num}, fieldname )
            focus_old = [getfield( data.CellA{cell_num}, fieldname),...
                focus];
        else
            focus_old = [focus];
        end
        
        data.CellA{cell_num} = setfield(data.CellA{cell_num}, fieldname, focus_old );
    end
    
    
    if disp_flag2
        figure(1);
        if cell_num
            plot( x,y, '.b' )
        else
            plot( x,y, '.r' )
            
        end
        
    end
    
end

%% Renormaze scores by cell
% Use the std of the model that include subtracking foci
for ii = 1:data.regs.num_regs
    
    if isfield( data.CellA{ii}, fieldname )
        
        xx = data.CellA{ii}.xx+crop;
        yy = data.CellA{ii}.yy+crop;
        
        im_ii = model_(yy,xx);
        
        mask_ii = data.CellA{ii}.mask;
        
        dI_ii = std( double(im_ii(mask_ii)));
        
        focus = getfield( data.CellA{ii}, fieldname  );
        num_loc = numel(focus);
        
        score_vec = [focus(:).score];
        [~,ord_sort] = sort( score_vec, 'descend' );
        focus = focus( ord_sort );
        
        for jj = 1:num_loc
            focus(jj).intensity_score = focus(jj).intensity_score*dI0/dI_ii;
            focus(jj).score           = focus(jj).score*dI0/dI_ii;
            
            focus(jj).shortaxis = ...
                (focus(jj).r-data.CellA{ii}.coord.rcm)*data.CellA{ii}.coord.e2;
            focus(jj).longaxis = ...
                (focus(jj).r-data.CellA{ii}.coord.rcm)*data.CellA{ii}.coord.e1;
        end
        
        
    else
        focus = focus0([]);
    end
    
    data.CellA{ii} = setfield( data.CellA{ii}, fieldname, focus );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = do_the_fit( ddd )
        
        gg = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4));
        
        %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        
        tmp = (double(imm)-gg);%.*mask_crop_for_fit(iy2,ix2);
        C = sum(tmp(logical(mkk)).^2);
        
        if disp_flag
            ddd(4);
            
            figure(3);
            clf;
            imshow( cat(3,ag([gg,imm]),ag([~mkk,~mkk]),0*[mkk,mkk]));
            '';
        end
    end

end






%%%%%%%%%%%%%%%%%%
%
% Gaussian Fit function

function C = g_fit_fun(X,Y,x0,y0,IG,I0,b1)

C = I0 + IG*exp( -((X-x0).^2+(Y-y0).^2)/(2*b1^2) );

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = makeImmortalRegions( data_mask, dirname_fluor1, CONST, disp_flag )

contents = dir( [dirname_fluor1,'*.tif'] );

numFluor = numel( contents );

s_im = size( data_mask.phase );

sum_im = zeros( s_im );

h = waitbar( 0, 'Making Sum image' );
for ii = 1:numFluor
    fluor = imread( [dirname_fluor1,contents(ii).name] );
    % compute the mean image of all the fluorescence frames
    sum_im = sum_im + double(fluor);
    
    waitbar( ii/numFluor, h );
end
close(h);

im0 = sum_im;

% Setup CONST calues for image processing 
CONST.SR.GausImgFilter_HighPass = fspecial('gaussian',141,10);
CONST.SR.GausImgFilter_LowPass3 = fspecial('gaussian',21,3);
CONST.SR.GausImgFilter_LowPass2 = fspecial('gaussian',21,2);
CONST.SR.GausImgFilter_LowPass1 = fspecial('gaussian',7,.75);

filt_disk =  strel( 'disk', CONST.count.intRad + 2);
diskInt   =  strel( 'disk', CONST.count.intRad );

% highpass image first to remove background laser intensity
im = double(im0) - imfilter(double(im0), CONST.SR.GausImgFilter_HighPass,'replicate');

% remove shot noise for dividing image into regions for fitting
im = imfilter(im,CONST.SR.GausImgFilter_LowPass1,'replicate');

% get mean/std intensity for the entire image
im_mean = mean(im(:));
im_std  = std(im(:));

% threshold the images at mean + std
im(im<(im_mean+2*im_std)) = im_mean+2*im_std;
im = im-(im_mean+im_std);

% make a mask of all regions above the cutoff.
mask = (im>0);
% remove stray pixels from the mask
mask_proc = bwmorph(bwmorph(mask,'erode'),'dilate');

% now divide regions that contain more than one fluorophore by
% performing a watershed.
ws = watershed(-im);
ws_tmp = and( data_mask.mask_cell,ws);
ws = bwlabel( ws_tmp );

% split locii
mask_proc2 = and( logical(ws), mask_proc);


regs_label_tmp = bwlabel( mask_proc2 );
num_regs = max( regs_label_tmp(:) );

%% make radius regions 
mask_new = false( s_im );
for ii = 1:num_regs

    mask_ii = (regs_label_tmp==ii);

    tmp = im  .*double(mask_ii);

    [junk, ind] = max(tmp(:));

    [y,x] = ind2sub(s_im,ind);
    
    mask_tmp = zeros( s_im );
    mask_tmp(y,x) = 1;
    
    mask_tmp = and( imdilate( mask_tmp, filt_disk ),...
        ws(y,x)==ws );
    
    
    
    mask_new = or( mask_new, mask_tmp);
    
end


% remove small regions
regs_label_tmp = bwlabel(mask_new);
num_regs = max( regs_label_tmp(:) );
props = regionprops( regs_label_tmp, {'Area','BoundingBox'} );

ind_tmp = find(and(and(...
drill( props, '.BoundingBox(3)' ) > 2,...
drill( props, '.BoundingBox(4)' ) > 2 ),...
[props.Area]' > CONST.findFocusSR.A_MIN ));

mask_new   = ismember( regs_label_tmp, ind_tmp);
mask_mod   = imerode( mask_new, diskInt );


regs_label = bwlabel(mask_new);
num_regs   = max( regs_label(:) );
props      = regionprops( regs_label, {'Area', 'BoundingBox'} );



data.sum_im      = sum_im;
data.sum_im_proc = im;
data.fluor_label = regs_label;

data.fluor_label_mod = double(mask_mod).*regs_label;

data.props       = props;
data.mask_mod    = mask_mod;



if disp_flag
    
    mask_cell = data_mask.mask_cell;
    
    outline_fr   = and(~mask_new,bwmorph(mask_new,'dilate'));
    outline_cell = and(~mask_cell,bwmorph(mask_cell,'dilate'));

    figure(1);
    clf;
    imshow( [cat(3,ag(sum_im),ag(sum_im),ag(sum_im)),...
             cat(3,ag(im)+0.3*ag(outline_fr),...
                   ag(im),...
                   ag(im)+0.6*ag(outline_cell))] );
end

end


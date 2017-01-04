%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = makeFluorRegions(data, CONST,dirname_xy, disp_flag );


%disp_flag = 1;


gaussR = CONST.trackLoci.gaussR;
crop   = CONST.trackLoci.crop;


tmp = struct(...
    'r', NaN, ...
    'score', NaN, ...
    'intensity_score', NaN, ...
    'intensity', NaN, ...
    'b', NaN );


persistent sdisk;
persistent gauss;
persistent Xpre;
persistent Ypre;
persistent crop_;



if isempty(sdisk)
    sdisk =  strel('disk',1);
    gauss  =  fspecial('gaussian',7,gaussR);
    %lap    = -fspecial('laplacian');
    %locusd =  imfilter(gauss, lap);
    
end


%%
%%---------------------------------------------------------%
% here we load the images and resize them
% the resizing is necessary because the fitting works better
% when there are many pixels per spot.
%
% im_mask_ is an enlarged cell mask to make the fitting
% more accurate when there are spots near the cell border
%---------------------------------------------------------


image = double(data.sum_im);


im_smooth = imfilter( image, gauss, 'replicate' );

im_mean = mean( im_smooth(:));
im_std = std( im_smooth(:));
im_thresh = im_mean + im_std;

im_smooth_th = im_smooth;

im_smooth_th( im_smooth<im_thresh) = im_thresh;
im_smooth_th = im_smooth_th-im_thresh;

L = watershed(-im_smooth_th);

mask_fluors = and(logical(L),data.mask_cell); 

fluor_label = bwlabel( mask_fluors );


data.im_smooth = im_smooth;

data.mask_fluors = mask_fluors;
data.fluor_label = fluor_label;

data.fluor_props = regionprops( fluor_label, {'Area', 'BoundingBox'} );



end





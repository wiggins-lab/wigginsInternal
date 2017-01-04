%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = makeImmortalRegions( data_mask, dir_cell, CONST, disp_flag )



num_dir_cell = numel( dir_cell );
im_sum_tot   = [];
im_sum_cell  = {};

for jj = 1:num_dir_cell
    contents = dir( [dir_cell{jj},'*.tif'] );
    
    if numel(contents) > 0
        numFluor = numel( contents );
        
        s_im = size( data_mask.phase );
        
        sum_im = zeros( s_im );
        
        h = waitbar( 0, 'Making Sum image' );
        
        
        for ii = 1:numFluor
            fluor = imread( [dir_cell{jj} ,contents(ii).name] );
            
            %fluor = fluor(200:300,200:300);
            % compute the mean image of all the fluorescence frames
            sum_im = sum_im + double(fluor);
            
            waitbar( ii/numFluor, h );
        end
        close(h);
        
        ranger1 = 148:150;
        ranger2 = 170:177;
        
        sum_im(ranger1,:) = median( sum_im(:) );
        sum_im(:,ranger2) = median( sum_im(:) );
        
        im_sum_max = max( sum_im(:) );
        im_sum_min = min( sum_im(:) );
        
        sum_im = (sum_im-im_sum_min)/(im_sum_max-im_sum_min);
        sum_im = sum_im/std(sum_im(:));
        
        im_sum_cell{jj} = sum_im;
        
        if isempty( im_sum_tot )
            im_sum_tot = sum_im;
        else
            im_sum_tot = im_sum_tot + sum_im;
        end
    end
end


sum_im = im_sum_tot;


im0 = sum_im; %(ranger1,ranger2);


%data_mask.mask_cell = data_mask.mask_cell(ranger1,ranger2);

% highpass image first to remove background laser intensity
CONST.SR.GausImgFilter_HighPass = fspecial('gaussian',141,10);
im = double(im0) - imfilter(double(im0), CONST.SR.GausImgFilter_HighPass,'replicate');


disp_flag_ = true;

if disp_flag_
    figure(1);
    clf;
    
   imshow( im, [] );
   hold on;
end

%% use our normal focus forming code to find loci and keep only foci that 
% look convincing

gf2 = fspecial( 'gaussian', 21, 2 );
opt =  optimset('MaxIter',25,'Display','off', 'TolX', 1/10);

locus = [];


%% make autofluor backers to use only if phase images exist of cells
if data_mask.phase_flag 
% calculate background of non cells
back_int = mean( double(im( ~data_mask.mask_cell)) );

% calculate the model for the background fluor.
backer =  imfilter( sqrt(double(bwdist( ~data_mask.mask_cell ))), gf2, 'replicate' );
end


% allow a max number of foci found per cell.
MAX_FOCI = 4;

% Calculate std in cells
Istd = std( double(im(logical( data_mask.mask_cell ))) );

% get size of image here.
sim = size( im0 );

disp_flag_ = true;
for ii = 1:data_mask.regs.num_regs
    
    [xx,yy]  = getBBpad( data_mask.regs.props(ii).BoundingBox, sim, 4 );
    
    mask_ii = data_mask.regs.regs_label(yy,xx) == ii;
    
    
    if data_mask.phase_flag
        locus_tmp = compSpotPosDevFmin4( im(yy,xx), im(yy,xx), backer(yy,xx), ...
        back_int, Istd, mask_ii, MAX_FOCI, CONST,1,disp_flag,2,opt);
    else
        locus_tmp = compSpotPosDevFmin2( im(yy,xx), im(yy,xx), mask_ii, MAX_FOCI, CONST,1,disp_flag,2,opt);  
    end
    
    if ~isempty( locus_tmp )
        
        
        % put locus positions in global coords
        num_spot = length(locus_tmp);
        for j = 1:num_spot;
            locus_tmp(j).r = locus_tmp(j).r + [xx(1)-1,yy(1)-1];   
            
            if disp_flag_
               text( locus_tmp(j).r(1)+3, locus_tmp(j).r(2), ...
                   num2str( locus_tmp(j).score,'%0.2f' ), 'Color', 'r' );
                
            end
        end
        
        
        % keep only foci which are within 25% of max score
        scores = [ locus_tmp(:).score];
        %keepers = scores > CONST.count.score_cut_rel*max(scores);
        keepers = scores > 2;
        locus_tmpk = locus_tmp(keepers);
        locus = [locus, locus_tmpk];
    end
end


%% Make a region for each locus using a watershed

% Setup CONST calues for image processing 
CONST.SR.GausImgFilter_LowPass3 = fspecial('gaussian',21,3);
CONST.SR.GausImgFilter_LowPass2 = fspecial('gaussian',21,2);
CONST.SR.GausImgFilter_LowPass1 = fspecial('gaussian',7,.75);

filt_disk =  strel( 'disk', CONST.count.intRad + 2);
diskInt   =  strel( 'disk', CONST.count.intRad );


% remove shot noise for dividing image into regions for fitting
im = imfilter(im,CONST.SR.GausImgFilter_LowPass1,'replicate');

% get mean/std intensity for the entire image
im_mean = mean(im(:));
im_std  = std(im(:));

% threshold the images at mean + std
%im(im<(im_mean+2*im_std)) = im_mean+2*im_std;
im = im-(im_mean);%+im_std);

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

% get keeper regions
keepers_regs = [];
for ii = 1:numel(locus)
    
    row = min( [max([1,round(locus(ii).r(2))]),sim(1)] );
    col = min( [max([1,round(locus(ii).r(1))]),sim(2)] );
    
    keep_tmp = regs_label_tmp(row,col);
    
    if ~keep_tmp 
        mask_new = false( sim );
        mask_new(row,col) = true;
        mask_new = bwdist( mask_new );
        mask_new( ~regs_label_tmp ) = nan;
        
        [~,ord] = min( mask_new(:) );
        keep_tmp = regs_label_tmp( ord );
    end

    
    keepers_regs = [keepers_regs,keep_tmp];
end
keepers_regs = keepers_regs(logical(keepers_regs));

regs_label_tmp = bwlabel( ismember(regs_label_tmp,keepers_regs) );
num_regs = max( regs_label_tmp(:) );

%% make radius regions 
mask_new = false( sim );
for ii = 1:num_regs

    mask_ii = (regs_label_tmp==ii);

    tmp = im  .*double(mask_ii);

    [~, ind] = max(tmp(:));

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
data.sum_im_cell = im_sum_cell;
data.sum_im_proc = im;
data.fluor_label = regs_label;

data.fluor_label_mod = double(mask_mod).*regs_label;

data.props       = props;
data.mask_mod    = mask_mod;

disp_flag = true;

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


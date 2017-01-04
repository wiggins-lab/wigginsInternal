%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = findFocusSR3(fluor, data_mask, data, CONST,opt,born,disp_flag)


persistent gauss;
% gaussian smoothing radius
gaussR       = CONST.findFocusSR.gaussR;
% gaussian smoothing matrix
if isempty(gauss)
    gauss  =  fspecial('gaussian',7,gaussR);
end


image      = double( fluor );
image_proc = imfilter( image, gauss, 'replicate' );
image_proc = image_proc - mean(image_proc(:));


% show fluor + phase composite image if disp_flag is true
if disp_flag
    figure(1);
    clf;
    imshow( cat(3,ag(~data_mask.mask_cell)*.0,ag(fluor),0.3*ag(data_mask.phase)) );
    
    figure(2)
    clf;
    imshow( fluor, [] );
end

% get the size of the phase image
sim = size( image );


num_regs = data_mask.regs.num_regs;
numFluor = data.numFluor;


%% init the trace variables if data has not been initialed
if ~isfield( data, 'regs' )
    % trace holds the traces
    
    % x is a vector of x positions
    trace0.x    = [];
    
    % y is a vector of y positions
    trace0.y    = [];
    
    % I is the amplitude of the fit gaussian
    trace0.I    = [];
    trace0.Isum    = [];
    trace0.IsumG    = [];
    trace0.IsumL    = [];
    
    
    % b is the standard deviation of the fit gaussian
    trace0.b    = [];
    
    
    % b is the standard deviation of the fit gaussian
    trace0.nn   = [];
    
    
    
    % born is the time this trace starts
    trace0.born =  -1;
    
    % n is the number of elements in the x,y,I,b vectors
    trace0.n    =  0;
    
    % n0 is the number of fits since the last fit with I
    % above the threshold CONST.findFocusSR.I_MIN
    trace0.n0   =  0;
    trace0.t0   =  4;
    trace0.Dt   =  1;
    trace0.I0G   =  [];
    trace0.I0L   =  [];
    
    % allow up to CONST.findFocusSR.MAX_TRACE_NUM traces
    trace(CONST.findFocusSR.MAX_TRACE_NUM) = trace0;
    for ii = 1:CONST.findFocusSR.MAX_TRACE_NUM
        % copy empty trace0 into the vector trace
        trace(ii) = trace0;
    end
    
    % flag showing if the trace is active or not.
    active = false( 1, CONST.findFocusSR.MAX_TRACE_NUM );
    
    % trace array
    regs0.trace    =  trace;
    
    % active array
    regs0.active   = active;
    
    % number of traces that have been fit
    regs0.numTrace =      0;
    
    % the offset of the cropped fluorescent image
    regs0.x0       =    NaN;
    regs0.y0       =    NaN;
    
    % the cell mask
    regs0.mask     =     [];
    
    % the vectors for the crop
    regs0.xx       =     [];
    regs0.yy       =     [];
    
    regs0.ndisk = [];

    % temporary
    regs0.xtmp = zeros( 1, CONST.findFocusSR.MAX_TRACE_NUM );
    regs0.ytmp = zeros( 1, CONST.findFocusSR.MAX_TRACE_NUM );
    
    regs0.ndisk = [];
    
    regs0.I0_mean = nan( 1, numFluor );
    regs0.I0_std  = nan( 1, numFluor );
    regs0.fluorProps =[];
    regs0.fluor_label = [];
    % Initialize a regs structure for each region
    regs(num_regs) = regs0;
    for ii = 1:num_regs
        
        % make the crop vectors with a pad of 4 pixels
        [xx,yy]  = getBBpad( data_mask.regs.props(ii).BoundingBox, sim, 4 );
        
        regs(ii) = regs0;
        
        % x0,y0 offset of the cropped image
        regs(ii).x0   = xx(1)-1;
        regs(ii).y0   = yy(1)-1;
        
        % store the crop vectors
        regs(ii).xx   = xx;
        regs(ii).yy   = yy;
        
        % make the mask for the ith region
        regs(ii).mask = (data_mask.regs.regs_label(yy,xx)==ii);
        
        
        
        tmp = data.fluor_label(yy,xx);
        tmp = unique(tmp( regs(ii).mask ));
        tmp = tmp(logical(tmp));
        tmp = reshape( tmp, [1,numel(tmp)] );
        
        regs(ii).fluorProps = data.props(tmp);
        
        regs(ii).fluor_label = data.fluor_label(yy,xx);
        
    end
    
    % copy the regs structure into the data structure
    data.regs = regs;
end

%% fit the focus positions, region by region.
for ii = 1:num_regs
%for ii = 2;
    
    mask_g = and( ~data.regs(ii).fluor_label, data.regs(ii).mask );
    xx = data.regs(ii).xx;
    yy = data.regs(ii).yy;

    im = fluor(yy,xx);

    data.regs(ii).I0_mean(born) = mean( double(im(mask_g)) );
    data.regs(ii).I0_std( born) = std(  double(im(mask_g)) );
    
    
    
    % get the fit for this frame, region ii
    data.regs(ii) = intCompSpotPos( image(yy,xx), image_proc(yy,xx), ...
        data.regs(ii), CONST, disp_flag, opt, born );
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function regs0 = intCompSpotPos(im_cell, im_proc, regs0, CONST, disp_flag, opt, born )


if ~exist('disp_flag', 'var');
    disp_flag = 0;
end


mask_cell = bwmorph(regs0.mask,'dilate' );;

channel = 1;

im_min = min(im_cell(:));
im_max = max(im_cell(:));



%
MAX_NUM_SPOT = CONST.findFocusSR.MAX_FOCUS_NUM;

% crop radius for the fit
crop         = CONST.findFocusSR.crop;

%
WS_CUT       = CONST.findFocusSR.WS_CUT;

% Number of amplitudes below threshold before
% a trace becomes inactive
MAX_OFF      = CONST.findFocusSR.MAX_OFF;

% Minimize intensity amplitude for the fit
I_MIN        = CONST.findFocusSR.I_MIN;

R_LINK       = CONST.findFocusSR.R_LINK;

A_MIN = 3;


% Make the variable only once
persistent sdisk;
persistent sdisk3;
persistent sdisk4;
persistent sdisk6;
persistent Xpre;
persistent Ypre;
persistent crop_;
persistent ndisk4;


if isempty(sdisk)
    % disk structure for mask dilation
    sdisk =  strel('disk',1);
    sdisk3 = strel( 'disk', 3);
    sdisk4 = strel( 'disk', 4);
    sdisk6 = strel( 'disk', 6);
    
    
    ndisk4 = sum(sum(sdisk4.getnhood));
    % crop vector for cropping images
    crop_ = -crop:crop;
    
    % 2D crop region array
    [Xpre, Ypre] = meshgrid(crop_,crop_);
end



if disp_flag
    disp( ['Time: ',num2str(born)] );
end

%%



% get image size
ss_cell = size(im_cell);

if disp_flag
    if ~exist('fig_num','var')
        fig_num = 7;
    end
end



% re-threshold the image to watershed it.
thresh  = WS_CUT;
im_thresh  = im_proc - thresh;
im_thresh(im_thresh<0) = 0;

% compute the watershed of the smoothed thresholded images
L      = watershed(-double(im_thresh));
L_mask = bwlabel( and( mask_cell, L));
wr     = (L==0);

num_fr = max(L_mask(:));


%%
% plots for debugging
if disp_flag
    figure(fig_num); clf;
    
    imshow([cat(3,autogain(im_cell),autogain(~mask_cell)*.2,autogain(wr)/2)],'InitialMagnification','fit');
end



%% Determine the fluor regions
% these are the positions of the maximum pixels in each
% region of the watershed
xn = zeros(1, num_fr);
yn = zeros(1, num_fr);

% this is the value of the max pixel
In = zeros(1, num_fr);

% this is the area of the watershed region
An = zeros(1, num_fr);

L_disk = zeros( ss_cell );
for j = 1:num_fr;
    mask_tmp = false( ss_cell );
    
    mask_j = (L_mask==j);
    
    [In(j),ind] = max(double(im_proc(:)).*double(mask_j(:)));
    [y,x] = ind2sub(ss_cell,ind);
    
    mask_tmp(y,x) = true;
    mask_tmp = and(mask_j, imdilate( mask_tmp, sdisk4));
    
    L_disk = L_disk + double(mask_tmp) * j;
    
    xn(j) = x;
    yn(j) = y;
    
    An(j) = sum( mask_tmp(:) );
end

L_disk = bwlabel(ismember( L_disk, find( and(An > A_MIN, In > I_MIN ))));
num_fr = max(L_disk(:));


actind = find( regs0.active );

% plots for debugging
if disp_flag
    figure(33); clf;
    imshow([cat(3,autogain(im_cell),autogain(~mask_cell)*.2,autogain(~L_disk)/2)],'InitialMagnification','fit');
   
    hold on;
    for ii = actind
       text( regs0.xtmp(ii),  regs0.ytmp(ii), num2str(ii), 'Color', 'y' );  
    end
end



%% compare to existing traces
actind = find( regs0.active );
nactind = numel( actind );

x_exist = round(regs0.xtmp(actind));
y_exist = round(regs0.ytmp(actind));

fr_ind = diag(L_disk( y_exist, x_exist ))';

% remove 2 to 1 mapping
[aa,bb1,cc] = unique( fr_ind, 'first' );
[aa,bb2,cc] = unique( fr_ind, 'last' );

ind_dup = find(and((bb1~=bb2),aa~=0));

for j = 1:numel(ind_dup)
   ind_tmp = actind(find( aa(ind_dup(j)) == fr_ind ));
 
   Ivals = drill( regs0.trace(ind_tmp),'.I(end)' );
   
   [junk, ind_tmp2] = max( Ivals );
   
   
  % 'hi'
   regs0.active(ind_tmp(ind_tmp~=ind_tmp(ind_tmp2))) = false;
end    

actind = find( regs0.active );
nactind = numel(actind);

x_exist = round(regs0.xtmp(actind));
y_exist = round(regs0.ytmp(actind));

fr_ind = diag(L_disk( y_exist, x_exist ))';

%%









%%

for j = 1:nactind
    
    % the trace isn't inside a fluor region, add one.
    if fr_ind(j) == 0
        
        mask_tmp = false( ss_cell );
        mask_tmp( y_exist(j), x_exist(j) ) = true;
        mask_tmp = imdilate( mask_tmp, sdisk4);
        mask_exist = bwmorph( logical(L_disk), 'dilate' );
        mask_tmp   = and( mask_tmp,and( mask_cell, ~mask_exist));
        
        num_reg_add = bwlabel( mask_tmp );
        num_reg_add = max( num_reg_add(:));
        
        [junk,ind] = max(double(im_proc(:)).*double(mask_tmp(:)));
        [y,x] = ind2sub(ss_cell,ind);
        
        if sum( mask_tmp(:) ) > A_MIN && num_reg_add==1;
            regs0.xtmp(actind(j)) = x;
            regs0.ytmp(actind(j)) = y;
            L_disk = L_disk - double(mask_tmp);
        else
            regs0.active(actind(j)) =false;
        end
        
    else
        ii = fr_ind(j);
        
        mask_tmp = (L_disk==ii);
        
        
        if sum( mask_tmp(:) ) < 3
            
            'hi'
        end
        
        [In(j),ind] = max(double(im_proc(:)).*double(mask_tmp(:)));
        [y,x] = ind2sub(ss_cell,ind);
        
        regs0.xtmp(actind(j)) = x;
        regs0.ytmp(actind(j)) = y;
        
    end
    
end





% If there aren't existing traces, add more.
fr_ind_all = 1:num_fr;
new_trace_fr_ind = fr_ind_all( ~ismember( fr_ind_all, fr_ind ));
n_new_trace_fr_ind = numel( new_trace_fr_ind );
numTrace = regs0.numTrace;
for j = 1:n_new_trace_fr_ind
    
    new_num = j+numTrace;
    regs0.active(new_num) = true;
    
    regs0.trace().born = born;
    
    ii = new_trace_fr_ind(j);
    mask_tmp = (L_disk==ii);
    [In(j),ind] = max(double(im_proc(:)).*double(mask_tmp(:)));
    [y,x] = ind2sub(ss_cell,ind);
    
    regs0.xtmp(new_num) = x;
    regs0.ytmp(new_num) = y;
    
    
end
numTrace = numTrace + n_new_trace_fr_ind;

regs0.numTrace = numTrace;

% re number regions after modifications.
L_disk = bwlabel(logical(L_disk));
props  = regionprops( L_disk, 'BoundingBox' );


%%
actind = find( regs0.active );
% plots for debugging
if disp_flag
    figure(44); clf;
    imshow([cat(3,autogain(im_cell),autogain(~mask_cell)*.2,autogain(~L_disk)/2)],'InitialMagnification','fit');
   
    hold on;
    for ii = actind
       text( regs0.xtmp(ii),  regs0.ytmp(ii), num2str(ii), 'Color', 'y' );  
    end
end


%% calculate the global backgorund value by removing foci
x000 = round(regs0.xtmp(actind));
y000 = round(regs0.ytmp(actind));

maskBack = false( ss_cell );
for ii = 1:numel( x000 )
    maskBack( y000(ii), x000(ii) ) = 1;
end

maskBack = and(~imdilate( maskBack, sdisk3 ),mask_cell);


I0G = mean( im_cell(maskBack) );





% Loop through the active regions
actind = find( regs0.active );
for ii = actind
    
    
    x00 = round(regs0.xtmp(ii));
    y00 = round(regs0.ytmp(ii));
    
    L_val = L_disk(y00,x00);
    
    [xx_fr,yy_fr] = getBB( props(L_val).BoundingBox );
    [XX_fr,YY_fr] = meshgrid(xx_fr,yy_fr);
    
    mask_fr = (L_disk(yy_fr,xx_fr)==L_val);
    im_fr   = im_cell(yy_fr,xx_fr);
    
    IG  = im_proc(y00,x00);
    b   = 2;
    
    ddd = [x00,y00,IG,b];
    
    
    xmin = xx_fr(1)-.5;
    xmax = xx_fr(end)+.5;
    
    ymin = yy_fr(1)-.5;
    ymax = yy_fr(end)+.5;
    
    dddMin = [xmin,ymin,im_min,1];
    dddMax = [xmax,ymax,im_max,4];
    
    [ddd] = lsqnonlin( @do_the_fit,ddd,dddMin,dddMax,opt);
    
    % extract the position, width, and intensity of the fitted gaussian
    
    if disp_flag
        im_T = g_fit_fun(XX_fr,YY_fr,ddd(1),ddd(2),ddd(3),I0G,ddd(4));
        
        figure(3)
        imshow( [im_fr,im_proc(yy_fr,xx_fr)+I0G,im_T], [] );
        colormap jet;
        drawnow;
    end
    
    
    x_ = ddd(1);
    y_ = ddd(2);
    I_ = ddd(3);
    b_ = ddd(4);
    
    x_r = round(x_);
    y_r = round(y_);
    
    % calculate the sum intensities
    mask_focus          = false( ss_cell );
    mask_focus(y_r,x_r) = true;
    mask_focus          = imdilate( mask_focus, sdisk4 );
    
    mask_local          = false( ss_cell );
    mask_local(y_r,x_r) = true;
    mask_local          = logical(double(imdilate( mask_local, sdisk6 ))...
        -double(mask_focus));
    
    I0L = mean( im_cell( mask_local ) );
    
    Isum     = sum( im_cell( mask_focus ) );
    IsumG    = sum( im_cell( mask_focus )-I0G );
    IsumL    = sum( im_cell( mask_focus )-I0L );
    
    
    % show the fit parameters
    if disp_flag
        disp(['Fit ',...
            ' x: ',num2str(round(ddd(1))),...
            ' y: ',num2str(round(ddd(2))),...
            ' I: ',num2str(round(ddd(3))),...
            ' b: ',num2str(ddd(4), '%1.2f')]);
    end
    
    % increment n for the fit
    regs0.trace(ii).n  = regs0.trace(ii).n  + 1;
    n_ = regs0.trace(ii).n;
    
    
    % copy (x,y,I, and b) into the next element in
    % the vector.
    regs0.trace(ii).x(n_)       =  x_ + regs0.x0;
    regs0.trace(ii).y(n_)       =  y_ + regs0.y0;
    regs0.trace(ii).I(n_)       =  I_;
    regs0.trace(ii).I0G(n_)       =  I0G;
    regs0.trace(ii).I0L(n_)       =  I0L;
    
    regs0.trace(ii).Isum(n_)    = Isum;
    regs0.trace(ii).IsumL(n_)   = IsumL;
    regs0.trace(ii).IsumG(n_)   = IsumG;
    regs0.trace(ii).b(n_)       =  b_;
    regs0.trace(ii).nn(n_)      =  born;
    regs0.ndisk                 =  ndisk4;
    
    % If the fit amplitude I is less than the cut off
    % I_MIN, increment n0
    if I_ < I_MIN
        regs0.trace(ii).n0 = regs0.trace(ii).n0 + 1;
        % set the (x,y) position to the previous (x,y)
        %         if numel(regs0.trace(ii).x) > 1
        %             regs0.trace(ii).x(end) = regs0.trace(ii).x((end-1));
        %             regs0.trace(ii).y(end) = regs0.trace(ii).y((end-1));
        %         end
        % but if I is above the cut off I_MIN, set n0 to 0
    else
        regs0.trace(ii).n0 = 0;
        regs0.xtmp(ii) = x_;
        regs0.ytmp(ii) = y_;
    end
    
    % if n0 > MAX_OFF, the have been too many fits that are
    % below threshold, assume bleached and make the trace
    % inactive
    if regs0.trace(ii).n0 > MAX_OFF
        regs0.active(ii) = false;
    end
    
    
    
    %% plots for debugging
    if disp_flag
        
        
        figure(2);
        hold on;
        plot(x_ + regs0.x0,y_ + regs0.y0,['r','.']);
        
        figure(1);
        hold on;
        plot(x_ + regs0.x0,y_  + regs0.y0,['r','.']);
        
        disp('');
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = do_the_fit( ddd )
        
        im_T = g_fit_fun(XX_fr,YY_fr,ddd(1),ddd(2),ddd(3),I0G,ddd(4));
        
        %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        
        tmp = (double(im_fr)-im_T).*mask_fr;
        C = tmp(:);
        
        %         figure(2);
        %         imshow( cat(3,ag([gg,imm]),ag([~mkk,~mkk]),0*[mkk,mkk]));
        %         '';
    end

end

%%%%%%%%%%%%%%%%%%
%
% Gaussian Fit function
%
%%%%%%%%%%%%%%%%%%

function C = g_fit_fun(X,Y,x0,y0,IG,I0,b1)

C = I0 + IG*exp( -((X-x0).^2+(Y-y0).^2)/(2*b1^2) );

end






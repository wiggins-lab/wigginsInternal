%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = findFocusSRc(fluor, data_mask, data, CONST,opt,born,disp_flag)



%% set up fluorescence

% mean fluor
fluorm  = mean( fluor(:) );

% standard deviation of fluor
fluors  = std(  double(fluor(:)) );

% threshold with mean + std 
%fluor( fluor< (fluorm+fluors) ) = fluorm+fluors;
%fluor = fluor - (fluorm+fluors);
%fluor( fluor< (fluorm) ) = fluorm;
fluor = double(fluor) - (fluorm);

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
sim = size( data_mask.phase );


num_regs = data_mask.regs.num_regs;


%% init the trace variables if data has not been initialed
if isempty( data )
    % trace holds the traces
    
    % x is a vector of x positions
    trace0.x    = [];
    
    % y is a vector of y positions
    trace0.y    = [];
    
    % I is the amplitude of the fit gaussian
    trace0.I    = [];
    trace0.Isum1    = [];
    trace0.Isum2    = [];
    trace0.Isum3    = [];

    trace0.Iint    = [];

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
    trace0.I0   =  0;
    
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
    
    % temporary 
    regs0.xtmp = zeros( 1, CONST.findFocusSR.MAX_TRACE_NUM );
    regs0.ytmp = zeros( 1, CONST.findFocusSR.MAX_TRACE_NUM ); 
    
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
    end
    
    % copy the regs structure into the data structure
    data.regs = regs;
end

%% fit the focus positions, region by region.
for ii = 1:num_regs
%for ii = 3;

    % im is the cropped fluorescence image.
    im   = fluor( data.regs(ii).yy, data.regs(ii).xx );

    % get the fit for this frame, region ii
    data.regs(ii) = intCompSpotPos( im, fluors, data.regs(ii), CONST, disp_flag, opt, born );
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function regs0 = intCompSpotPos(image, ims, regs0, CONST, disp_flag, opt, born )


if ~exist('disp_flag', 'var');
    disp_flag = 0;
end


im_mask = bwmorph(regs0.mask,'dilate' );;
channel = 1;


% 
MAX_NUM_SPOT = CONST.findFocusSR.MAX_FOCUS_NUM;

% gaussian smoothing radius
gaussR       = CONST.findFocusSR.gaussR;

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



% Make the variable only once
persistent sdisk;
persistent sdisk3;
persistent sdisk4;
persistent sdisk6;
persistent gauss;
persistent Xpre;
persistent Ypre;
persistent crop_;

if isempty(sdisk)
    % disk structure for mask dilation 
    sdisk =  strel('disk',1);
    sdisk3 = strel( 'disk', 3);
    sdisk4 = strel( 'disk', 4); 
    sdisk6 = strel( 'disk', 6); 

    % gaussian smoothing matrix
    gauss  =  fspecial('gaussian',7,gaussR);
    
    % crop vector for cropping images
    crop_ = -crop:crop;
    
    % 2D crop region array
    [Xpre, Ypre] = meshgrid(crop_,crop_);  
end



if disp_flag
    disp( ['Time: ',num2str(born)] );
end

%%


image    = double(image);
im_mask  = logical(im_mask);

im_mask_dilated = imdilate(im_mask, sdisk);
im_mask_dilated(1,:) = 0;
im_mask_dilated(end,:) = 0;
im_mask_dilated(:,1) = 0;
im_mask_dilated(:,end) = 0;


% get image size
imsize = size(image);

if disp_flag
    if ~exist('fig_num','var')
        fig_num = 7;
    end
end

%---------------------------------------------------------
%
% pad the images to facilitate cropping around each spot
%
% there are two masks: one for watersheding and one for fitting.
% the fitting mask is slightly wider to accomodate spots near
% the cell boundary.
%
%---------------------------------------------------------

im = zeros( imsize(1)+2*crop,imsize(2)+2*crop);
mk = logical(im);
mk_for_fit = logical(im);

mk(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = im_mask;
%mk(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = im_mask_dilated;
im(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = image;
mk_for_fit(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = im_mask_dilated;


% mean and variance of the intensity within the cell
dI = ims;

% here we segment the fluorescent image to obtain regions
% around each spot.



% add pad to the outside of the image
im_loc = zeros( imsize(1)+2*crop,imsize(2)+2*crop);
im_loc(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = ...
    imfilter( (image), gauss, 'replicate' );

% re-threshold the image to watershed it.
thresh  = WS_CUT;
im_loc  = im_loc - thresh;
im_loc(im_loc<0) = 0;

% compute the watershed of the smoothed thresholded images
L = watershed(-double(im_loc));

wr = (L == 0);
num_seg = max((L(:)));

% these are the positions of the maximum pixels in each 
% region of the watershed
xn = zeros(num_seg,1);
yn = zeros(num_seg,1);

% this is the value of the max pixel
fn = zeros(num_seg,1);

% this is the area of the watershed region
an = zeros(num_seg,1);

% zero out pixels of the pad region
LL = L;
LL(~mk_for_fit) = 0;

%%
% plots for debugging
if disp_flag
    figure(fig_num); clf;
    if channel == 1
        % imshow([cat(3,autogain(im),autogain(~mk)*.2,autogain(wr)/2),...
        %         cat(3,autogain(im),autogain(~mk)*0,autogain(im_loc))],'InitialMagnification','fit')
        imshow([cat(3,autogain(im),autogain(~mk)*.2,autogain(wr)/2)],'InitialMagnification','fit');
    elseif channel == 2
        imshow(cat(3,autogain(wr)/2,0.5*autogain(im.*mk),autogain(im.*mk)),'InitialMagnification','fit')
    else
        imshow(cat(3,autogain(im.*mk)*.75,autogain(im.*mk)*.75,autogain(wr)/2),'InitialMagnification','fit')
    end
    %    drawnow;
    hold on
    %   axis equal
end


%%
% get maxima from each region
if disp_flag
        figure(8);
        clf;
        imshow( im_loc, [] );
end

for j = 1:num_seg;
    
    mask = (LL==j);
    
    [maxValue,y] = max(double(im_loc(:)).*double(mask(:)));
    [y,x] = ind2sub(size(im_loc),y);
    
    xn(j) = x;
    yn(j) = y;
    %
    %     bb=size(image);
    %     tmp = im(max(y-crop,1):min(y+crop,bb(1)),max(x-crop,1):min(x+crop,bb(2)));
    fn(j)=maxValue;
    
    an(j) = sum(mask(:));
    
    if disp_flag

        hold on;
        plot( x, y, 'r.' );
        disp('');
    end
end

%% take the pixels above the cut
AREA_MIN = 4;

indcut = and( (fn > I_MIN), (an > AREA_MIN) );
fn = fn(indcut)';
yn = yn(indcut)';
xn = xn(indcut)';
an = an(indcut)';


num_cut = numel(yn);


%% display the guess using the max pixel as a guess
if disp_flag
    for gg = 1:numel(yn)
        disp(['Pre ',...
            ' x: ',num2str(round(xn(gg))),...
            ' y: ',num2str(round(yn(gg))),...
            ' I: ',num2str(round(fn(gg)))]);
    end
end


% if born == 26
%     'hi'
% end


%% find x that aren't there and then add them to the list

actind = find( regs0.active );
nactind = numel( actind );

[ind1_no_match, ind1_match, indA_match] = ...
    SR2( xn - crop + regs0.x0, yn - crop + regs0.y0, drill(regs0.trace(actind),'.x(end)'), ...
    drill(regs0.trace(actind),'.y(end)'), CONST.findFocusSR.R_LINK);

n_ind1_no_match = numel( ind1_no_match );
n_ind1_match    = numel( ind1_match );

if n_ind1_no_match
    new_ind = 1:n_ind1_no_match;
    for kkk = new_ind
        regs0.trace(kkk+regs0.numTrace).born = born;
        %regs0.trace(kkk+regs0.numTrace).x    =    xn(mind(kkk)) - crop + regs0.x0;
        %regs0.trace(kkk+regs0.numTrace).y    =    yn(mind(kkk)) - crop + regs0.y0;
        regs0.trace(kkk+regs0.numTrace).I     =    fn(ind1_no_match(kkk));
        regs0.trace(kkk+regs0.numTrace).b     =    2;
        regs0.active(kkk+regs0.numTrace)      =    true;
        regs0.trace(kkk+regs0.numTrace).Isum1 =    0*fn(ind1_no_match(kkk));
        regs0.trace(kkk+regs0.numTrace).Isum2 =    0*fn(ind1_no_match(kkk));
        regs0.trace(kkk+regs0.numTrace).Isum3 =    0*fn(ind1_no_match(kkk));
         
        regs0.xtmp(kkk+regs0.numTrace)    =    xn(ind1_no_match(kkk)) - crop + regs0.x0;
        regs0.ytmp(kkk+regs0.numTrace)    =    yn(ind1_no_match(kkk)) - crop + regs0.y0;
        regs0.trace(kkk+regs0.numTrace).I0     =   0*fn(ind1_no_match(kkk));
    end
    
    actind = [actind, new_ind+regs0.numTrace];
    regs0.numTrace = regs0.numTrace + n_ind1_no_match;
end

if n_ind1_match
    try
        regs0.xtmp( actind(indA_match) ) = xn(ind1_match) - crop + regs0.x0;
        regs0.ytmp( actind(indA_match) ) = yn(ind1_match) - crop + regs0.y0;
    catch
        'hi'
    end
end

%%
%background intensity set to zero
%I0 = sum( double(im(:)).*double(mk(:)) )/sum( double(mk(:)) );
%I0 = 0;


x000 = round(regs0.xtmp(actind) + crop - regs0.x0);
y000 = round(regs0.ytmp(actind) + crop - regs0.y0);

maskBack = 0*mk_for_fit;

for ii = 1:numel( x000 )
    maskBack( y000(ii), x000(ii) ) = 1; 
end

maskBack = and(~imdilate( maskBack, sdisk3 ),mk_for_fit);

%imshow( maskBack, [] );

I0 = mean( im(maskBack) );
%I0 = ims;

for ii = actind
    
    %x = round(regs0.trace(ii).x(end) + crop - regs0.x0);
    %y = round(regs0.trace(ii).y(end) + crop - regs0.y0);
    x00 = round(regs0.xtmp(ii) + crop - regs0.x0);
    y00 = round(regs0.ytmp(ii) + crop - regs0.y0);
    
    
    %[regs0.trace(ii).x(:) + crop - regs0.x0,regs0.trace(ii).y(:) + crop - regs0.y0]
    
    mask_crop_for_fit = and(bwmorph(LL==LL(y00,x00),'dilate'),mk_for_fit);
    
    % this is the initial guess for the fit parameters
    
    %[size(im),y,x,regs0.y0,regs0.x0]
    IG  = im(y00,x00);
    b   = 2;
    
    ddd = [x00,y00,IG,b];
    
    
    x = min([imsize(2),max([1,x00])]);
    y = min([imsize(1),max([1,y00])]);
    
   
    
    y = min([imsize(1),max([1,y00])]);
    
    
    % lower bound
    ix2 = x+crop_;
    iy2 = y+crop_;

    xx2 = Xpre(logical(iy2),logical(ix2))+x;
    yy2 = Ypre(logical(iy2),logical(ix2))+y;
        
    ix2 = ix2(logical(ix2));
    iy2 = iy2(logical(iy2));
    
    % here we actually do the fit
    
    mkk = mask_crop_for_fit(iy2,ix2);
    imm = im(iy2,ix2);
     
    xmin = ix2(1,1);
    xmax = ix2(end,end);  
    
    ymin = iy2(1,1);
    ymax = iy2(end,end); 
    
    
    dddMin = [xmin,ymin,-max(image(:))*1.2,1];
    dddMax = [xmax,ymax,max(image(:))*1.2,4];
    
    [ddd] = lsqnonlin( @do_the_fit,ddd,dddMin,dddMax,opt);
    
    % extract the position, width, and intensity of the fitted gaussian
   
   if disp_flag
        immT = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4));
       
        figure(3)
        imshow( [imm,immT], [] );
        colormap jet;
        drawnow;
    end
    
    
    x_ = ddd(1);
    y_ = ddd(2);
    I_ = ddd(3);
    b_ = ddd(4);
    
    
    % calculate the sum intensities
    mask_sum4      = false( size( im ) ); 
    mask_sum4(y,x) = true;
    mask_sum4      = imdilate( mask_sum4, sdisk4 );

    mask_sum6      = false( size( im ) ); 
    mask_sum6(y,x) = true;
    mask_sum6      = logical(double(imdilate( mask_sum6, sdisk6 ))...
        -double(mask_sum4));

    Isum2    = sum( im( mask_sum4 )-I0 );
    Isum3    = sum( im( mask_sum4 )...
        -mean( im( mask_sum6 ) ));

    
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
    
    
    
    
    x_other = regs0.xtmp(actind) + crop - regs0.x0;
    y_other = regs0.ytmp(actind) + crop - regs0.y0;

    dr = sqrt( (x_-x_other).^2 + (y_-y_other).^2 );
    
    [junk, ind_min_dr ] = min( dr );
    
    
    % if the fit (x_,y_) is crop away from the
    % previous fit, set (x,y) to the last fit 
    % position and set the amplitude I = 0.
    ddr = R_LINK;
    
    flag_reset = (abs(x_-x) > ddr) || ...
       (abs(y_-y) > ddr) || ...
       (actind(ind_min_dr)~=ii);
   
    if flag_reset
   
        x_ = x;
        y_ = y;
        I_ = im(y,x)-I0;
        %I_ = 0;
        b  = 1.5;
    end
    
    % copy (x,y,I, and b) into the next element in 
    % the vector.
    regs0.trace(ii).x(n_)       =  x_ - crop + regs0.x0;
    regs0.trace(ii).y(n_)       =  y_ - crop + regs0.y0;
    regs0.trace(ii).I(n_)       =  I_;
    regs0.trace(ii).I0(n_)       =  I0;
    regs0.trace(ii).Isum1(n_)   = pi*2*b_^2*I_;
    regs0.trace(ii).Isum2(n_)   = Isum2;
    regs0.trace(ii).Isum3(n_)   = Isum3;
    regs0.trace(ii).b(n_)       =  b_;
    regs0.trace(ii).nn(n_)      =  born;

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
    end
    
    % if n0 > MAX_OFF, the have been too many fits that are
    % below threshold, assume bleached and make the trace 
    % inactive
    if regs0.trace(ii).n0 > MAX_OFF
        regs0.active(ii) = false;
    end
   
    regs0.xtmp(ii) = regs0.trace(ii).x(end);
    regs0.ytmp(ii) = regs0.trace(ii).y(end);
    
    %% plots for debugging
    if disp_flag
        
        cc = 'w';
        
        figure(fig_num);
        
        if ~flag_reset
        plot(x_,y_,[cc,'.']);        
        text(x_+2,y_,['n: ', num2str(ii),' I: ',num2str(I_,3)],'Color',cc);
        end
        
        
        plot(x,y,'go');
        
        
        figure(2);
        hold on;
        plot(x_ - crop + regs0.x0,y_ - crop + regs0.y0,['r','.']);
        
        figure(1);
        hold on;
        plot(x_ - crop + regs0.x0,y_ - crop + regs0.y0,['r','.']);
        
        disp('');
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = do_the_fit( ddd )
        
        gg = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4));
        
        %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        
        tmp = (double(imm)-gg).*mask_crop_for_fit(iy2,ix2);
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






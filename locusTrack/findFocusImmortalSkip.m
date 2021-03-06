%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = findFocusImmortal(fluor, data_mask, data, CONST, opt, born, disp_flag)

persistent gauss;

% gaussian smoothing radius
gaussR       = CONST.findFocusSR.gaussR;



if isempty(gauss)
    
    % gaussian smoothing matrix
    gauss  =  fspecial('gaussian',7,gaussR);
    
end



%% set up fluorescence
% add pad to the outside of the image
fluor_proc = imfilter( fluor, gauss, 'replicate' );


% mean fluor
fluorm  = mean( fluor_proc(:) );

% standard deviation of fluor
fluors  = std(  double(fluor_proc(:)) );


numFluorS = data.numFluorS;
numFluor  = data.numFluor;


% threshold with mean + std
%fluor( fluor< (fluorm+fluors) ) = fluorm+fluors;
%fluor = fluor - (fluorm+fluors);
%fluor( fluor< (fluorm) ) = fluorm;


%fluor_proc = double(fluor_proc) - (fluorm);







% show fluor + phase composite image if disp_flag is true
if disp_flag
    figure(1);
    clf;
    imshow( cat(3,ag(~data_mask.mask_cell)*.0,ag(fluor),0.3*ag(data_mask.phase)) );
    
    figure(2)
    clf;
    imshow( fluor_proc, [] );
end

% get the size of the phase image
sim = size( data_mask.phase );


num_regs = data_mask.regs.num_regs;


%% init the trace variables if data has not been initialed
if ~isfield( data, 'regs' )
    % trace holds the traces
    
    % x is a vector of x positions
    trace0.x    = nan( 1, numFluorS );
    trace0.y    = nan( 1, numFluorS );
    trace0.x_s  = nan( 1, numFluorS );
    trace0.y_s  = nan( 1, numFluorS );
    
    % I is the amplitude of the fit gaussian
    trace0.I    = nan( 1, numFluorS );
    trace0.Isum    = nan( 1, numFluorS );
    trace0.IsumG    = nan( 1, numFluorS );
    trace0.IsumL    = nan( 1, numFluorS );
    
    trace0.Iint    = nan( 1, numFluorS );
    
    % b is the standard deviation of the fit gaussian
    trace0.b    = nan( 1, numFluorS );
    
    
    % b is the standard deviation of the fit gaussian
    trace0.nn   = nan( 1, numFluorS );
    
    
    
    % born is the time this trace starts
    trace0.born =  -1;
    
    % n is the number of elements in the x,y,I,b vectors
    trace0.n    =  0;
    
    % n0 is the number of fits since the last fit with I
    % above the threshold CONST.findFocusSR.I_MIN
    trace0.n0   =  0;
    trace0.t0   =  4;
    trace0.Dt   =  1;
    trace0.I0   =  nan( 1, numFluorS );
    
    trace0.mask   =  [];
    trace0.xx     =  [];
    trace0.yy     =  [];
    
    trace0.active   =  false( 1, numFluorS );
    
    
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
    regs0.xtmp = zeros( 1, numFluorS );
    regs0.ytmp = zeros( 1, numFluorS );
    regs0.mask =[];
    regs0.fluorNum = [];
          
    regs0.fluorProps = [];
    regs0.fluor_label = [];
    
    regs0.ndisk = [];
    
    regs0.I0_mean = [];
    regs0.I0_std  = [];
    
    regs0.I0_meanS = nan( 1, numFluorS );
    regs0.I0_stdS  = nan( 1, numFluorS );
    
    regs0.traceS = trace0 ;
    
    regs0.kymo = [];
    regs0.kymo_mask = [];
    regs0.kymo_ss = [];
    
    regs0.im_cell = cell( [1, numFluor] );

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
        
        regs(ii).numTrace = numel(tmp);
        regs(ii).fluorNum = tmp;
        
        if  regs(ii).numTrace 
        regs(ii).traceS( regs(ii).numTrace ) = trace0;
        end
        
        regs(ii).fluorProps = data.props(tmp);
        
        regs(ii).fluor_label = data.fluor_label(yy,xx);
        
        % initialize kymo 
        regs(ii).kymo_mask =  imrotate( regs(ii).mask, -data_mask.regs.props(ii).Orientation, 'bilinear'  );
        regs(ii).kymo_ss = size( regs(ii).kymo_mask );
        ss_tmp = size( regs(ii).kymo_mask );
        regs(ii).kymo      = zeros( ss_tmp(2), numFluorS );
        
        for jj = 1:regs(ii).numTrace
            kk = tmp( jj );
            
            [xx,yy]  = getBBpad( data.props(kk).BoundingBox, sim, 2 );
            regs(ii).traceS(jj) = trace0;
            regs(ii).traceS(jj).mask = (data.fluor_label(yy,xx)==kk);
            
            regs(ii).traceS(jj).xx = xx;
            regs(ii).traceS(jj).yy = yy;
        end
    end
    
    % copy the regs structure into the data structure
    data.regs = regs;
end



%% fit the focus positions, region by region.
for ii = 1:num_regs

    
    mask_g = and( ~data.regs(ii).fluor_label, data.regs(ii).mask );
    xx = data.regs(ii).xx;
    yy = data.regs(ii).yy;

    im = fluor(yy,xx);

    data.regs(ii).I0_meanS(born) = mean( double(im(mask_g)) );
    data.regs(ii).I0_stdS( born) = std(  double(im(mask_g)) );

    
    % get the fit for this frame, region ii
    data.regs(ii) = intCompSpotPos( fluor, fluor_proc, data.regs(ii), CONST, disp_flag, opt, born );
    
   %% Take care of kymo
 kymo_tmp = imrotate( fluor_proc(yy,xx), -data_mask.regs.props(ii).Orientation, 'bilinear'  );
 
 I0g = data.regs(ii).I0_meanS(born );
 kymo_tmp(kymo_tmp<I0g) = I0g;
 kymo_tmp = kymo_tmp - I0g;
 
 kymo_tmp = kymo_tmp.*double(data.regs(ii).kymo_mask);
 data.regs(ii).kymo(:,born) = sum( kymo_tmp );
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function regs0 = intCompSpotPos(image, imageS, regs0, CONST, disp_flag, opt, n_ )

% mask for global intnensity
xx  = regs0.xx;
yy  = regs0.yy;





%%
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
persistent ndisk4;

if isempty(sdisk4)
    % disk structure for mask dilation

    sdisk4 = strel( 'disk', 4);
    
    ndisk4 = sum(sum(sdisk4.getnhood));
end



if disp_flag
    disp( ['Time: ',num2str(n_)] );
end

%%

numTrace =regs0.numTrace;



I0g  = regs0.I0_meanS(n_);
I0gs = regs0.I0_stdS(n_);



 
 


 
 
%% Do fitting here

for ii = 1:numTrace;
    
    mask_ii    = regs0.traceS(ii).mask;
    
    
    
    xx = regs0.traceS(ii).xx;
    yy = regs0.traceS(ii).yy;
    
    im_ii      = double(image( yy,xx));
    
    im_proc = imageS(yy,xx);
    s_im_proc = size( im_proc );
    
    [maxValue,ind] = max(double(im_proc(:)).*double(mask_ii(:)));
    [y,x]          = ind2sub( s_im_proc, ind);
    
    
    x0 = x+xx(1)-1;
    y0 = y+yy(1)-1;
    
    
    %[size(im),y,x,regs0.y0,regs0.x0]
    IG  = abs(double(imageS(y0,x0))-I0g);
    IG = max( [IG, 1] );
    
    
    b   = 2;
    
    ddd = double([x0,y0,IG,b]);
    
    
    [xx2,yy2] = meshgrid( xx,yy );
    
    
    
    dddMin = double([xx(1),   yy(1),   0, 1]);
    dddMax = double([xx(end), yy(end),  1e6,     4]);
    
    gg = [];
    
    
    if ~isnan( I0g )
        [ddd] = lsqnonlin( @do_the_fit,ddd,dddMin,dddMax,opt);
    end
    % extract the position, width, and intensity of the fitted gaussian
    
    if disp_flag
        immT = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0g,ddd(4));
        
        figure(3)
        imshow( [im_ii,im_proc,immT], [] );
        colormap jet;
        drawnow;
    end
    
    
    x_ = ddd(1);
    y_ = ddd(2);
    I_ = ddd(3);
    b_ = ddd(4);
    
    
    
    
    x_int = round( x_ );
    y_int = round( y_ );
    
    
    
    if ~mask_ii(y_int-yy(1)+1,x_int-xx(1)+1) || ...
            any( ddd(1:2)==dddMax(1:2) ) || ...
            any( ddd(1:2)==dddMin(1:2) )
        
        if n_ > 1
            x_ = regs0.traceS(ii).x(n_-1);
            y_ = regs0.traceS(ii).y(n_-1);
            
        else
            x_ = nan;
            y_ = nan;
        end
        
        I_ = 0;
        b_ = nan;
        
    end
    
    x_int = round( x_ );
    y_int = round( y_ );
    
    if isnan(x_int) || isnan(y_int)
        x_int = x0;
        y_int = y0;
    end
    
    
    mask_focus = false( s_im_proc );
    mask_focus(y_int-yy(1)+1,x_int-xx(1)+1) = true;
    
    mask_focus = imdilate( mask_focus, sdisk4);
    
    mask_local = and( mask_ii,~mask_focus);
    
    I0l = mean( im_ii( mask_local ) );
    
    Isum  = sum( im_ii( mask_focus ) ); 
    IsumG = sum( im_ii( mask_focus )-I0g );
    
    IsumL = sum( im_ii( mask_focus )-I0l );
    
    if isnan( IsumL )
        IsumL = IsumG;
    end
    
    % show the fit parameters
    if disp_flag
        disp(['Fit ',...
            ' x: ',num2str(round(ddd(1))),...
            ' y: ',num2str(round(ddd(2))),...
            ' I: ',num2str(round(ddd(3))),...
            ' b: ',num2str(ddd(4), '%1.2f')]);
    end
    
    
    regs0.traceS(ii).x(n_)       =  x_;
    regs0.traceS(ii).y(n_)       =  y_;
    regs0.traceS(ii).I(n_)       =  I_;
    regs0.traceS(ii).I0(n_)      =  I0g;
    regs0.traceS(ii).Isum(n_)    = Isum;
    regs0.traceS(ii).IsumG(n_)   = IsumG;
    regs0.traceS(ii).IsumL(n_)   = IsumL;
    regs0.traceS(ii).b(n_)       =  b_;
    regs0.traceS(ii).nn(n_)      =  n_;
    regs0.traceS(ii).n           =  n_;
    regs0.ndisk                 =  ndisk4;
        
    regs0.traceS(ii).active(n_)  = ( IsumG > 2*I0gs*sqrt(ndisk4));
     
    if  regs0.traceS(ii).active(n_) || n_ == 1
        regs0.traceS(ii).x_s(n_) = regs0.traceS(ii).x(n_);
        regs0.traceS(ii).y_s(n_) = regs0.traceS(ii).y(n_);
    else
        regs0.traceS(ii).x_s(n_) = regs0.traceS(ii).x_s(n_-1);
        regs0.traceS(ii).y_s(n_) = regs0.traceS(ii).y_s(n_-1);
        regs0.traceS(ii).x(n_)   = regs0.traceS(ii).x(n_-1); 
        regs0.traceS(ii).y(n_)   = regs0.traceS(ii).y(n_-1);
    end
    
    
    %% plots for debugging
    if disp_flag
        
        cc = 'w';
        
        
        figure(2);
        hold on;
        plot(x_,y_,['r','.']);
        
        figure(1);
        hold on;
        plot(x_,y_ ,['r','.']);
        
        disp('');
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = do_the_fit( ddd )
        
        gg = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0g,ddd(4));
        
        %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        
        tmp = (double(im_ii)-gg).*mask_ii;
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






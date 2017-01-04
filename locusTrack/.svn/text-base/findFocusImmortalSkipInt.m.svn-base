%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = findFocusImmortalSkipInt(fluor, data_mask, data, CONST, opt, born, n_s, disp_flag)

persistent gauss;
persistent sdisk4;

fluor_m = medfilt2( fluor, [3,3], 'symmetric' );

% gaussian smoothing radius
gaussR       = CONST.findFocusSR.gaussR;



if isempty(gauss)
    
    % gaussian smoothing matrix
    gauss  =  fspecial('gaussian',7,gaussR);
    sdisk4 = strel( 'disk', CONST.count.intRad);
end



%% set up fluorescence
% add pad to the outside of the image
fluor_proc = imfilter( fluor, gauss, 'replicate' );


% mean fluor
fluorm  = mean( fluor_proc(:) );

% standard deviation of fluor
fluors  = std(  double(fluor_proc(:)) );


numFluor = data.numFluor;


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
if born == 1
    % trace holds the traces
    
    % x is a vector of x positions
    trace0.x    = nan( 1, numFluor );
    trace0.x_s  = nan( 1, numFluor );
    % y is a vector of y positions
    trace0.y    = nan( 1, numFluor );
    trace0.y_s  = nan( 1, numFluor );
    % I is the amplitude of the fit gaussian
    trace0.I        = nan( 1, numFluor );
    trace0.Isum     = nan( 1, numFluor );
    trace0.IsumG    = nan( 1, numFluor );
    trace0.IsumL    = nan( 1, numFluor );
    
    trace0.Iint     = nan( 1, numFluor );
    
    % I is the amplitude of the fit gaussian
    trace0.Im         = nan( 1, numFluor );
    trace0.Im_sum     = nan( 1, numFluor );
    trace0.Im_sumG    = nan( 1, numFluor );
    trace0.Im_sumL    = nan( 1, numFluor );
    
    trace0.Im_int     = nan( 1, numFluor );
    
    
    % b is the standard deviation of the fit gaussian
    trace0.b    = nan( 1, numFluor );
    
    
    % b is the standard deviation of the fit gaussian
    trace0.nn   = nan( 1, numFluor );
    
    
    
    % born is the time this trace starts
    trace0.born =  -1;
    
    % n is the number of elements in the x,y,I,b vectors
    trace0.n    =  0;
    
    % n0 is the number of fits since the last fit with I
    % above the threshold CONST.findFocusSR.I_MIN
    trace0.n0   =  0;
    trace0.t0   =  4;
    trace0.Dt   =  1;
    trace0.I0   =  nan( 1, numFluor );
    
    trace0.mask   =  [];
    trace0.xx     =  [];
    trace0.yy     =  [];
    
    trace0.active =  false( 1, numFluor );
    
    trace0.micro_im  = cell( [1, numFluor] );
    trace0.micro_imm = cell( [1, numFluor] );
    
    for ii = 1:num_regs
        
        % make the crop vectors with a pad of 4 pixels
        
        if   data.regs(ii).numTrace      
                data.regs(ii).trace( data.regs(ii).numTrace ) = trace0;
        end
        
        
        for jj = 1:data.regs(ii).numTrace
            
            
            data.regs(ii).trace(jj).mask = data.regs(ii).traceS(jj).mask;
            
            data.regs(ii).trace(jj).xx   = data.regs(ii).traceS(jj).xx;
            data.regs(ii).trace(jj).yy   = data.regs(ii).traceS(jj).yy;
        end
    end
    
    % copy the regs structure into the data structure
end



%% fit the focus positions, region by region.
for ii = 1:num_regs
    
    xx = data.regs(ii).xx;
    yy = data.regs(ii).yy;
    
    % make mask for background subtraction
    mask_g = false( size(fluor) );
    numTrace = data.regs(ii).numTrace;
    for jj = 1:numTrace;
        
        x_ = round(data.regs(ii).traceS(jj).x_s(n_s));
        y_ = round(data.regs(ii).traceS(jj).y_s(n_s));
        if ~isnan( y_ ) && ~isnan( x_ )
            mask_g(y_,x_) = true;
        end
    end
    mask_g = imdilate( mask_g, sdisk4 );
    mask_g = and( ~mask_g(yy,xx), data.regs(ii).mask );
    
    %mask_g = and( ~data.regs(ii).fluor_label, data.regs(ii).mask );
    
    im  = fluor(yy,xx);
    imm = fluor_m(yy,xx);

    data.regs(ii).I0_mean(born) = mean( double(im(mask_g)) );
    data.regs(ii).I0_std( born) = std(  double(im(mask_g)) );
    
    data.regs(ii).I0m_mean(born) = mean( double(imm(mask_g)) );
    data.regs(ii).I0m_std( born) = std(  double(imm(mask_g)) );
    
    % get the fit for this frame, region ii
    data.regs(ii) = intCompSpotPos( fluor, fluor_m, fluor_proc, data.regs(ii), ...
        CONST, disp_flag, opt, born, n_s );
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function regs0 = intCompSpotPos(image, image_m, imageS, regs0, CONST, disp_flag,...
    opt, n_, n_s )


if ~exist('disp_flag', 'var');
    disp_flag = 0;
end


im_mask = bwmorph(regs0.mask,'dilate' );;
channel = 1;






% Make the variable only once
persistent sdisk4;
persistent ndisk4;
persistent fdisk4;


if isempty(sdisk4)
    % disk structure for mask dilation
    sdisk4 = strel(    'disk', CONST.count.intRad);
    fdisk4 = fspecial( 'disk', CONST.count.intRad);

    ndisk4 = sum(sum(sdisk4.getnhood));
end



if disp_flag
    disp( ['Time: ',num2str(n_)] );
end



%%

numTrace =regs0.numTrace;

% mask for global intnensity
xx  = regs0.xx;
yy  = regs0.yy;

I0g  = regs0.I0_mean(n_);
Im0g = regs0.I0m_mean(n_);

%disp(I0g)

for ii = 1:numTrace;
    
    mask_ii    = regs0.trace(ii).mask;
    
    
    
    xx = regs0.trace(ii).xx;
    yy = regs0.trace(ii).yy;
    
    im_ii      = double(image( yy,xx));
    imm_ii     = double(image_m( yy,xx));

    s_im_proc = size( im_ii );
    
    
    x_ = regs0.traceS(ii).x_s(n_s);
    y_ = regs0.traceS(ii).y_s(n_s);
    
    
    x_int = round( x_ );
    y_int = round( y_ );
    
    if isnan(x_int)  || isnan(y_int)
        IsumG = nan;
        IsumL = nan;
        Isum  = nan;
        
        ImsumG = nan;
        ImsumL = nan;
        Imsum  = nan;
        
        ndisk = nan; 
    else
        mask_focus = false( size(image) );
        mask_ii_   = mask_focus;
        
        mask_focus(y_int,x_int) = true;
        
        mask_focus = imdilate( mask_focus, sdisk4);
        
        mask_ii_(yy,xx) = mask_ii;
        
        mask_local = and( mask_ii_,~mask_focus);
        
        % from raw im
        I0l = mean( image( mask_local ) );
        
        Isum  = sum( image( mask_focus ) );
        IsumG = sum( image( mask_focus )-I0g );
        
        IsumL = sum( image( mask_focus )-I0l );
        
        if isnan( IsumL )
            IsumL = IsumG;
        end
        
        % from med image
        Im0l = mean( image_m( mask_local ) );
        
        Imsum  = sum( image_m( mask_focus ) );
        ImsumG = sum( image_m( mask_focus )-Im0g );
        
        ImsumL = sum( image_m( mask_focus )-Im0l );
        
        
        ndisk = sum(mask_focus(:));
    end
    
    regs0.trace(ii).x(n_)       = regs0.traceS(ii).x(n_s);
    regs0.trace(ii).y(n_)       =  regs0.traceS(ii).y(n_s);
    
    regs0.trace(ii).x_s(n_)     = x_;
    regs0.trace(ii).y_s(n_)     = y_;
    
    regs0.trace(ii).I(n_)       =  regs0.traceS(ii).I(n_s);
    
    regs0.trace(ii).I0(n_)      = I0g;
    regs0.trace(ii).Isum(n_)    = Isum;
    regs0.trace(ii).IsumG(n_)   = IsumG;
    regs0.trace(ii).IsumL(n_)   = IsumL;
    
    regs0.trace(ii).Im0(n_)      = Im0g;
    regs0.trace(ii).Imsum(n_)    = Imsum;
    regs0.trace(ii).ImsumG(n_)   = ImsumG;
    regs0.trace(ii).ImsumL(n_)   = ImsumL;
    
    regs0.trace(ii).b(n_)       =  regs0.traceS(ii).b(n_s);
    regs0.trace(ii).nn(n_)      =  n_;
    regs0.trace(ii).n           =  n_;
    regs0.trace(ii).active(n_)  =  regs0.traceS(ii).active(n_s);
    regs0.ndisk(ii)             =  ndisk;
    
    regs0.trace(ii).micro_im{n_}  =  im_ii;
    regs0.trace(ii).micro_imm{n_} =  imm_ii;

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
end

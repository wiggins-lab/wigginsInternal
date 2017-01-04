%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function data = findFocusImmortalSkipIntW(fluor, data_mask, data, CONST, opt, t, n_s, disp_flag)

persistent gauss;
persistent sdisk4;
persistent fdisk4;

fluor_m = medfilt2( fluor, [3,3], 'symmetric' );

% gaussian smoothing radius
gaussR       = CONST.findFocusSR.gaussR;
intRad       = CONST.count.intRad;


if isempty(gauss)
    
    % gaussian smoothing matrix
    gauss  =  fspecial( 'gaussian', 7,gaussR);
    sdisk4 =  strel(    'disk',     intRad);
    fdisk4 =  fspecial( 'disk',     intRad)*pi*intRad^2;
end


npixels = sum( fdisk4(:) );

%% set up fluorescence
% add pad to the outside of the image
fluor_proc = imfilter( fluor, gauss,  'replicate' );
fluor_disk = imfilter( fluor, fdisk4, 'replicate' );


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

ss = size( fluor_disk );





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
if t == 1
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
    
    for ii = 1:num_regs
        
        % make the crop vectors with a pad of 4 pixels
        
        for jj = 1:data.regs(ii).numTrace
            data.regs(ii).trace(jj) = trace0;
        end
        
        
        for jj = 1:data.regs(ii).numTrace
            
            
            data.regs(ii).trace(jj).mask = data.regs(ii).traceS(jj).mask;
            
            data.regs(ii).trace(jj).xx   = data.regs(ii).traceS(jj).xx;
            data.regs(ii).trace(jj).yy   = data.regs(ii).traceS(jj).yy;
        end
    end
    
    % copy the regs structure into the data structure
end

%%
disp_flag = false;

if disp_flag
    figure(11);
    clf;
    
    imshow( fluor_disk, [] );
    hold on;
    colormap jet;
    
    
    figure(1);
    clf;
    
    imshow( fluor, [] );
    hold on;
    colormap jet;
end

%% fit the focus positions, region by region.
for ii = 1:num_regs
    
    xx = data.regs(ii).xx;
    yy = data.regs(ii).yy;
    
    
    % Compute integrated fluor
    mask_g = false( size(fluor) );
    numTrace = data.regs(ii).numTrace;
    
    data.regs(ii).im_cell{t} = fluor(yy,xx);
    
    for jj = 1:numTrace;
        
        mask_jj = (data.fluor_label(yy,xx) == data.regs(ii).fluorNum(jj));
        
        if any( mask_jj(:) )
            
            fluor_dt = fluor_disk(yy,xx);
            fluor_dt( ~mask_jj ) = nan;
            
            
            if CONST.count.moveFlag
                [Isum, ord_jj] = max( fluor_dt(:) );
                [y,x]          = ind2sub( size(fluor_dt), ord_jj);
                y_ = y + yy(1)-1;
                x_ = x + xx(1)-1;
            else
                y_   = max( [1,min( [round( data.regs(ii).traceS(jj).y(1) ), ss(1)] )] );
                x_   = max( [1,min( [round( data.regs(ii).traceS(jj).x(1) ), ss(2)] )] );
                Isum = fluor_disk( y_,x_ );  
                y    = y_ - yy(1)+1;
                x    = x_ - xx(1)+1;
            end
            
            data.regs(ii).trace(jj).Isum(t) = Isum;
            data.regs(ii).trace(jj).x(t)    = x_;
            data.regs(ii).trace(jj).y(t)    = y_;
            
            mask_g(y_,x_) = true;
            
            if disp_flag
                figure(11)
                plot( x_, y_, '.k' );
                figure(1)
                plot( x_, y_, '.k' );
            end
        end
        
    end
end

mask_g_ = imdilate( mask_g, sdisk4 );

for ii = 1:num_regs
    
    xx = data.regs(ii).xx;
    yy = data.regs(ii).yy;
    
    fluor_ii = fluor( yy, xx );
    % Compute integrated fluor
    numTrace = data.regs(ii).numTrace;
    
    mask_ii = and( data.regs(ii).mask, ~mask_g_(yy,xx) );
    
    Isum0G  = mean( double(fluor_ii(mask_ii) ) );
    dIsum0G = std(  double(fluor_ii(mask_ii) ) );

    
    data.regs(ii).I0_mean(t) = mean( double(fluor_ii(mask_ii) ) );
    data.regs(ii).I0_std( t) = std(  double(fluor_ii(mask_ii) ) );
    
    
    for jj = 1:numTrace;
        mask_jj = (data.fluor_label(yy,xx) == data.regs(ii).fluorNum(jj));
        Isum0L = sum( double(fluor_ii(mask_jj) ) )/sum( double(mask_jj(:)) ) ;
        
        data.regs(ii).trace(jj).IsumG(t) = data.regs(ii).trace(jj).Isum(t) - ...
            Isum0G*npixels;
        
        data.regs(ii).trace(jj).IsumL(t) = data.regs(ii).trace(jj).Isum(t) - ...
            Isum0L*npixels;
        
        data.regs(ii).trace(jj).I0(t)     = Isum0G;
        data.regs(ii).trace(jj).dI0(t)    = dIsum0G;
        
        xx_jj = data.regs(ii).trace(jj).xx;
        yy_jj = data.regs(ii).trace(jj).yy;
        
        data.regs(ii).trace(jj).micro_im{t}  =  fluor( yy_jj, xx_jj );
        
        
        data.regs(ii).trace(jj).x_s(t)     = data.regs(ii).trace(jj).x(t);
        data.regs(ii).trace(jj).y_s(t)     = data.regs(ii).trace(jj).y(t);
        
        data.regs(ii).trace(jj).I(t)       =  data.regs(ii).traceS(jj).I(n_s);
        
        data.regs(ii).trace(jj).b(t)       =  data.regs(ii).traceS(jj).b(n_s);
        data.regs(ii).trace(jj).nn(t)      =  t;
        data.regs(ii).trace(jj).n           =  t;
        data.regs(ii).trace(jj).active(t)  =  true;
        data.regs(ii).ndisk(jj)    =  npixels;
        
    end
end




end


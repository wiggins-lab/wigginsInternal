%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [data, with_intensity, framesWithIntensity] = ...
    findFocusT6(framesWithIntensity, with_intensity, fluor, data, CONST, opt, disp_flag, time)




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


%% set up fluorescence

fluor_smooth = imfilter( fluor, gauss, 'replicate' );

mask_bg  = data.mask_bg;

%filter noise
filt = fspecial('laplacian');
toRemove = imfilter(fluor_smooth, filt, 'replicate');
fluor_smooth = fluor_smooth - toRemove;

pix = fluor_smooth(mask_bg);
pix(pix == 0) = [];


% mean fluor
fluorm  = mean( fluor_smooth(:) );

% standard deviation of fluor
fluors  = std(  double(fluor_smooth(:)) );

% mean fluor subtracting minimum pixel
fluorm_cell = mean(pix);

% standard deviation of fluor subtracting minimum pixel
fluors_cell = std(double(pix));



% threshold with mean + std
%fluor( fluor< (fluorm+fluors) ) = fluorm+fluors;
%fluor = fluor - (fluorm+fluors);
%fluor( fluor< (fluorm) ) = fluorm;
fluor = double(fluor) - (fluorm);

% show fluor + phase composite image if disp_flag is true
if disp_flag
    figure(1);
    clf;
    imshow( cat(3,ag(~data.mask_cell)*.0,ag(fluor),0.3*ag(data.phase)) );
    
    figure(2)
    clf;
    imshow( fluor, [] );
end

% get the size of the phase image
sim = size( data.phase );


num_regs = data.regs.num_regs;


%% init the trace variables if data has not been initialed
%if ~isfield( data, 'cell' )
    % trace holds the traces
    
    
    % x is a vector of x positions
    trace0.x    = nan(1,data.num_t);
    
    % y is a vector of y positions
    trace0.y    = nan(1,data.num_t);
    
    % I is the amplitude of the fit gaussian
    trace0.I    = nan(1,data.num_t);
    trace0.Isum1    = nan(1,data.num_t);
    trace0.Isum2    = nan(1,data.num_t);
    trace0.Isum3    = nan(1,data.num_t);
    
    trace0.Iint    = nan(1,data.num_t);
    
    % b is the standard deviation of the fit gaussian
    trace0.b    = nan(1,data.num_t);
    
    trace0.I0   =  nan(1,data.num_t);
    
    trace0.fl_reg_num = nan;
    
    trace0.mask = [];
    
    trace0.x0 = [];
    trace0.y0 = [];
    
    trace0.xx = [];
    trace0.yy = [];
    
    trace0.x_offset = [];
    trace0.y_offset = [];
    
    trace0.sim = [];
    
    trace0.maxValue = [];
    
    
    for ii = 1:data.regs.num_regs
        
        tmp = unique(data.fluor_label(data.regs.regs_label(:)==ii));
        
        tmp = tmp( logical( tmp ));
        
        data.cell(ii).fluor_reg_num = tmp;
        data.cell(ii).reg_num = numel(tmp);
        
        
        for jj = 1:data.cell(ii).reg_num
            
            if time == 1
                data.cell(ii).trace(jj) = trace0;
            end
            
            
            kk = data.cell(ii).fluor_reg_num(jj);
            
            % make the crop vectors with a pad of 4 pixels
            [xx,yy]  = getBBpad( data.fluor_props(kk).BoundingBox, sim, 4 );
            
            data.cell(ii).trace(jj).xx = xx;
            data.cell(ii).trace(jj).yy = yy;
            
            data.cell(ii).trace(jj).x_offset = xx(1)-1;
            data.cell(ii).trace(jj).y_offset = yy(1)-1;
            
            
            mask = ...
                bwmorph((data.fluor_label(yy,xx)==kk),'dilate' );
            
            
            data.cell(ii).trace(jj).mask = mask;
            
            data.cell(ii).trace(jj).sim = size( mask );
            
            im_smooth = fluor_smooth(yy,xx);
            [maxValue,y] = max(double(im_smooth(:)).*double(mask(:)));
            
            [y0,x0] = ind2sub(size(im_smooth),y);
            
            data.cell(ii).trace(jj).y0 = y0 + data.cell(ii).trace(jj).y_offset;
            data.cell(ii).trace(jj).x0 = x0 + data.cell(ii).trace(jj).x_offset;
            
            
            
        end
        
    end
%end



%% loop through all regions corresponding to cells

for kk = 1:num_regs
    
    % in each cell loop through all regions
    
    for jj = 1:data.cell(kk).reg_num
        
        
        
        
        xx = data.cell(kk).trace(jj).xx;
        yy = data.cell(kk).trace(jj).yy;
        
        mask = data.cell(kk).trace(jj).mask;
        
        
        im_smooth = fluor_smooth(yy,xx);
        s_im_proc = size( im_smooth );
        
        [maxValue,ind] = max(double(im_smooth(:)).*double(mask(:)));
        [y00,x00]          = ind2sub( s_im_proc, ind);
        
        
        image = fluor( yy, xx);
        
        image    = double(image);
        
        image_mean = mean( image(:) );
        
        im = image;
        im( im< image_mean) =image_mean;
        im = im - image_mean;
        
        x_offset = xx(1)-1;
        y_offset = yy(1)-1;
        
        x0 = x00 + x_offset;
        y0 = y00 + y_offset;
        
        
        %x0 = data.cell(kk).trace(jj).x0;
        %y0 = data.cell(kk).trace(jj).y0;
        
        %x00 = data.cell(kk).trace(jj).x0 - x_offset;
        %y00 = data.cell(kk).trace(jj).y0 - y_offset;
        
        IG  = im(y00,x00);
        b   = 2;
        
        ddd = [x0,y0,IG,b];
        
        
        % local coords for crop image
        ix2 = x00+crop_;
        iy2 = y00+crop_;
        
        sim_loc = data.cell(kk).trace(jj).sim;
        
        ix2(ix2>sim_loc(2)) = 0;
        ix2(ix2<1)          = 0;
        iy2(iy2>sim_loc(1)) = 0;
        iy2(iy2<1)          = 0;
        
        
        xx2 = Xpre(logical(iy2),logical(ix2))+x0;
        yy2 = Ypre(logical(iy2),logical(ix2))+y0;
        
        ix2 = ix2(logical(ix2));
        iy2 = iy2(logical(iy2));
        
        % here we actually do the fit
        
        mkk = mask(iy2,ix2);
        imm = im(iy2,ix2);
        
        
        
        
        xmin = ix2(1,1)+x_offset;
        xmax = ix2(end,end)+x_offset;
        
        ymin = iy2(1,1)+y_offset;
        ymax = iy2(end,end)+y_offset;
        
        
        
        
        dddMin = [xmin,ymin,-max(image(:))*1.2,1];
        dddMax = [xmax,ymax,max(image(:))*1.2,4];
        
        I0 = 0;
        
        % set threshold, so if maxValue is too small, don't fit (save
        % time!)
        if maxValue > (fluorm_cell + 3*fluors_cell)
            %disp( 'do it' );
            
            with_intensity(kk) = kk; %if this is interesting cell, mark it
            framesWithIntensity(kk,time) = kk;
            
            [ddd] = lsqnonlin( @do_the_fit,ddd,dddMin,dddMax,opt);
            
            % extract the position, width, and intensity of the fitted gaussian
            
            %             if disp_flag
            %                 immT = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4));
            %
            %                 figure(3)
            %                 imshow( [imm,immT], [] );
            %                 colormap jet;
            %                 drawnow;
            %             end
            
            
            x_ = ddd(1);
            y_ = ddd(2);
            I_ = ddd(3);
            b_ = ddd(4);
            
            
            % show the fit parameters
            %         if disp_flag
            %             disp(['Fit ',...
            %                 ' x: ',num2str(round(ddd(1))),...
            %                 ' y: ',num2str(round(ddd(2))),...
            %                 ' I: ',num2str(round(ddd(3))),...
            %                 ' b: ',num2str(ddd(4), '%1.2f')]);
            %         end
            
            % increment n for the fit
            
            
            
        else
            %disp( 'skip' );
            
            
            x_ = nan;
            y_ = nan;
            I_ = nan;
            b_ = nan;
            
            I0 = nan;
        end
        
        data.cell(kk).trace(jj).x(    time)       =  x_;
        data.cell(kk).trace(jj).y(    time)       =  y_;
        data.cell(kk).trace(jj).I(time)           =  I_;
%         
%         if kk == 12
%             disp(jj)
%             data.cell(kk).trace(jj).I
%         end
        
        data.cell(kk).trace(jj).I0(time)          =  I0;
        
        
        data.cell(kk).trace(jj).b(time)       =  b_;
        
        %% plots for debugging
        if disp_flag
            
            cc = 'w';
            
            %             figure(2);
            %
            % %            imshow( fluor, [] );
            %             hold on;
            %
            %                 plot(x_,y_,[cc,'.']);
            %                 text(x_+2,y_,['n: ', num2str(ii),' I: ',num2str(I_,3)],'Color',cc);
            %
            %
            %             plot(x,y,'go');
            
            
            figure(2);
            hold on;
            plot(x_ ,y_,['r','.']);
            
            
        end
        
    end
        
end

% tot = 0;
% for iii = 1:data.cell(12).reg_num
%     tot = tot + nansum(data.cell(12).trace(iii).I);
% end
% tot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = do_the_fit( ddd )
        
        gg = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4));
        
        %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        
        tmp = (double(imm)-gg).*double(mkk);
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






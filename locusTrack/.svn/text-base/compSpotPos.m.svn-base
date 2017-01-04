function [spot] = compSpotPos(image, im_mask, crop, channel, disp_flag, fig_num)

spot = [];

if ~exist('fig_num')
    fig_num = 7;
end

%---------------------------------------------------------
%
% here we load the images and resize them
% the resizing is necessary because the fitting works better
% when there are many pixels per spot.
%
% im_mask_ is an enlarged cell mask to make the fitting
% more accurate when there are spots near the cell border
%---------------------------------------------------------


image = double(image);
im_mask = double(im_mask);

im_mask_dilated = imdilate(im_mask, strel('disk',1));
im_mask_dilated(1,:) = 0;
im_mask_dilated(end,:) = 0;
im_mask_dilated(:,1) = 0;
im_mask_dilated(:,end) = 0;

% image = imresize(image,2,'nearest');
% im_mask = imresize(im_mask,2,'nearest');
% im_mask_dilated = imresize(im_mask_dilated,2,'nearest');

imsize = size(image);

if ~exist('disp_flag');
    disp_flag = 0;
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
mk = im;
mk_for_fit = im;

mk(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = im_mask;
im(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = image;
mk_for_fit(crop+1:crop+imsize(1),crop+1:crop+imsize(2)) = im_mask_dilated;


%---------------------------------------------------------
%
% mean and variance of the intensity within the cell
%
%---------------------------------------------------------

Ibar = mean(im(logical(mk)));
dI = std(im(logical(mk)));

imsize = size(im);


%---------------------------------------------------------
%
% here we segment the fluorescent image to obtain regions
% around each spot.
%
%---------------------------------------------------------



% im_smoothed = imfilter(im_smoothed,fspecial('gaussian',5,1/2),'replicate');
im_smoothed = imfilter(im,fspecial('gaussian',7,1),'replicate');
im_smoothed = double(im_smoothed).*mk;

L = watershed(-double(im_smoothed));

wr = (L == 0);
num_seg = max((L(:)));

xn = zeros(num_seg,1);
yn = zeros(num_seg,1);
fn = zeros(num_seg,1);
an = zeros(num_seg,1);

spot_num = 0;

LL = L.*mk_for_fit;

%---------------------------------------------------------
%
% here we look at each watershed region and record
% the intensity of its spot.
%
%---------------------------------------------------------

for j = 1:num_seg;
    
    mask = (LL==j).*mk;
    im_smoothed_tmp = double(mask).*im_smoothed;
    
    [maxValue,y] = max(im_smoothed_tmp(:));
    [y,x] = ind2sub(size(im),y);
    
    xn(j) = x;
    yn(j) = y;
    %
    %     bb=size(image);
    %     tmp = im(max(y-crop,1):min(y+crop,bb(1)),max(x-crop,1):min(x+crop,bb(2)));
    fn(j)=maxValue;
    
    an(j) = sum(mask(:));
    
end

%---------------------------------------------------------
%
% plots for debugging
%
%---------------------------------------------------------
if disp_flag
    figure(fig_num); clf;
    if channel == 1
        imshow(cat(3,autogain(im),autogain(~mk)*.2,autogain(wr)/2),'InitialMagnification','fit')
    elseif channel == 2
        imshow(cat(3,autogain(wr)/2,0.5*autogain(im.*mk),autogain(im.*mk)),'InitialMagnification','fit')
    else
        imshow(cat(3,autogain(im.*mk)*.75,autogain(im.*mk)*.75,autogain(wr)/2),'InitialMagnification','fit')
    end
%    drawnow;
    hold on
 %   axis equal
end


%---------------------------------------------------------
%
% now we fit a gaussian to locate the spot within each region
%
%---------------------------------------------------------


%background intensity set to zero
%I0 = sum( double(im(:)).*double(mk(:)) )/sum( double(mk(:)) );
I0 = Ibar;

% thresholds to determine which spots are real
fnn = (fn - Ibar)/dI;

%%%%%%%%%%%%%%%%%%%%%%

% MODIFICATION: ONLY FIT THE BRIGHTEST TWO WATERSHEDS

%%%%%%%%%%%%%%%%%%%%%%%%

MIN_AREA = crop^2;



[fnn_, ord] = sort(fnn, 'descend');
an = an(ord);

iis = find(  (fnn_>0) .* (an > MIN_AREA)  );



iis = iis(1:min([2,numel(iis)]));
iis = reshape(iis,1,numel(iis));

if isempty(iis)
    return;
end


for ii = iis
    
    i = ord(ii);
    spot_crop = 2;
    mk(yn(i)+(-spot_crop:spot_crop), xn(i)+(-spot_crop:spot_crop)) = 0;
    
end


Ibar = mean(im(~~mk));
dI = std(im(~~mk));


if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clf;
% back = autogain(im_smoothed);
% wsag = autogain(mask_crop_for_fit);
% imshow( cat(3,autogain(im_smoothed),autogain(mask_crop_for_fit),autogain(mask_crop_for_fit)), [] );
% hold on;
% for ii = iis'
%     i = ord(ii);
%     plot(xn(i),yn(i),'wo');
%     
% end
% 'hi'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end







opt =  optimset('MaxIter',10,'Display','off', 'TolX', 1/10,'TolFun', 1e-8);


for ii = iis
    
    i = ord(ii);
    
    x = xn(i);
    y = yn(i);
    
    %     mask_crop = double(L==i).*mk;
    mask_crop_for_fit = double(LL==i);
    
    
    
    % this is the initial guess for the fit parameters
    IG = im(y,x)-I0;
    b = 4;
    bmin = 1;
    bmax = 6;
    
    ddd = [x,y,IG,b,100,b];
    
    % lower bound
    lddd = [x-crop,y-crop,0,bmin,-100,bmin];
    
    % upper bound
    uddd = [x+crop,y+crop,100000,bmax,100,bmax];
    
    % make the coordinates for the cropped image
    ix2 = x+(-crop:crop);
    iy2 = y+(-crop:crop);
    [xx2,yy2] = meshgrid(ix2,iy2);
    
    
    %---------------------------------------------------------
    %
    % here we actually do the fit
    %
    %---------------------------------------------------------
    

    mkk = mask_crop_for_fit(iy2,ix2);
    imm = im(iy2,ix2);
    Ibar_local = min(imm(:));
    
    IbarMax = max([Ibar_local,Ibar]);
    %imm = uint8(imm-IbarMax)+IbarMax;
    %imm = uint8(imm-Ibar);
    

    
    [ddd, resn, res, exit_flag, output, lambda, jac] = lsqnonlin( @do_the_fit, ddd, lddd, uddd, opt );
    
    
    %          if findstr(output.algorithm, 'line')
    % ''
    %          end
    
    %---------------------------------------------------------
    %
    % extract the position, width, and intensity of the fitted gaussian
    %
    %---------------------------------------------------------
    
    xpos = ddd(1);
    ypos = ddd(2);
    
    
    sz = size(L);
    
    position_score = LL(min(max(1,floor(ypos)),sz(1)), min(max(1,floor(xpos)),sz(2)));
    position_score = (position_score==i);
    

    
    b1 = ddd(4);
    b2 = ddd(6);
    IG = ddd(3);
    theta = ddd(5);
    
    B = max(ddd([4,6]));
    

    Iamp = ddd(3);
   
    Ifit = g_fit_fun(xx2,yy2,xpos,ypos,IG,I0,b1,theta,b2).*mkk;
    Itotal = sum(double((Ifit(:)-Ibar(:)).*mkk(:)));
    Ivar = dI*sqrt(sum(double(logical(mkk(:)))));
    
    if disp_flag
    figure(2);
    imshow( autogain(([Ifit, imm.*mkk])) ,'InitialMagnification','fit');
    '';
    end
    
    score_i = Itotal./Ivar-1;
    score = position_score*score_i;
    
    
    
    
    intensity_score = position_score*((ddd(3))/dI -1);
    raw_intensity = ddd(3);
    
    
   
    
    %---------------------------------------------------------
    %
    % only record the spot if the fit looks good
    %
    %---------------------------------------------------------
    
    if score > 0
        
        spot_num = spot_num+1;
        
        rr(:,i) = [xpos,ypos];
        spot(spot_num).r = [xpos-crop,ypos-crop];
        spot(spot_num).score = score;
        spot(spot_num).intensity_score = intensity_score;
        spot(spot_num).intensity = raw_intensity;
        spot(spot_num).b = B;
        score_array(spot_num) = score;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % factor of 1/2 here because we resized the images
        %
        
        %             spot(spot_num).r = spot(spot_num).r/2;
        %             spot(spot_num).b = spot(spot_num).b/2;
        
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    %---------------------------------------------------------
    %
    % plots for debugging
    %
    %---------------------------------------------------------
    
    if disp_flag
        
        if score > 0;
            cc = 'w';
        else
            cc = 'b';
        end
        
        figure(fig_num);
        plot(xpos,ypos,[cc,'.']);
        text(xpos+2,ypos,['s:', num2str(score,3)],'Color',cc);
        %text(xpos+2,ypos,[num2str(round(intensity_score*100)/100),'::', num2str(round(score*100)/100)],'Color',cc);
        plot(xpos+b1*[0,sin(theta)],ypos+b1*[0,cos(theta)],[cc,':']);
        plot(xpos+b2*[0,cos(theta)],ypos+b2*[0,-sin(theta)],[cc,':']);
        
        
        %             n=0:.01:1;
        %             theta = n*2*pi;
        %             yt = sin(theta);
        %             xt = cos(theta);
        %             plot(xpos+B*yt,ypos+B*xt,[cc,':']);
    end
    
end


% Sort spots by score so that first spot has the highest score.

if exist('score_array')
    [a,b] = sort(score_array,'descend');
    spot = spot(b);
end



    function C = do_the_fit( ddd )
        
        gg = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4),ddd(5),ddd(6));
        
        %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
        
        tmp = (double(imm)-gg);%.*mask_crop_for_fit(iy2,ix2);
        C = tmp(logical(mkk));
        
    end

end

%%%%%%%%%%%%%%%%%%
%
% Gaussian Fit function
%
%%%%%%%%%%%%%%%%%%

function C = g_fit_fun(X,Y,x0,y0,IG,I0,b1,theta,b2)

rot = [[cos(theta), sin(theta)];[-sin(theta), cos(theta)]];
KK = rot*diag([b1^-2,b2^-2])*rot'/2;
K11 = KK(1,1);
K12 = KK(1,2);
K22 = KK(2,2);

C = I0 + IG*exp( -(K11*(X-x0).^2+2*K12*(Y-y0).*(X-x0)+K22*(Y-y0).^2) );

end






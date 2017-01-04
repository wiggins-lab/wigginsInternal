%% Model of cytoplasmic fluor in cell. Fit cell by cell.
% numc is channel number
% function adds the field cytoX, where cyto is the global cytofluor model 
% for channel X = numc.
function [ data ] = intFindFociCyto( data, CONST, numc )



debug_flag = true;


disp_flag = false;
disp_flag2 = false;

fieldname = ['locus',num2str(numc)];

% Get images out of the structures.
image0 = double(getfield( data, ['fluor',num2str(numc)] ));

% Do background subtract
hg = fspecial( 'gaussian' ,210,30 );
image = image0 - imfilter( image0, hg, 'replicate' );
image_ = image;

%image = image0;
%image = medfilt2( image0, [3,3], 'symmetric' );


dI0 = std( image(data.mask_bg) ); 


image_med = medfilt2( image, [3,3], 'symmetric' );

if debug_flag
    figure(1);
    clf;
    imshow( image_med, [] );
    hold on;
end



hg = fspecial( 'gaussian' ,14,1 );
images = imfilter( image_med, hg, 'replicate' );

rad = 3;
A   = pi*rad^2;
hd = fspecial( 'disk' , rad)*A;
imaged = imfilter( image, hd, 'replicate' );


mask_mod = bwmorph( data.mask_bg, 'dilate', 1 );

L = watershed( -images );
Lp = logical(double(L).*double(mask_mod));

fr_label = bwlabel( Lp );
props    = regionprops( fr_label, {'BoundingBox'} );
num_regs = numel( props );

imsize = size( image );


focus0.r               = [nan,nan];
focus0.score           = nan;
focus0.intensity_score = nan;
focus0.intensity       = nan;
focus0.b               = 1;
focus0.error           = nan;
focus0.shortaxis       = nan;
focus0.longaxis        = nan;

%% Loop through the regions and fit loci in each one

max_med_ii  = nan( [1,num_regs] );
max_disk_ii = nan( [1,num_regs] );
x_ii        = nan( [1,num_regs] );
y_ii        = nan( [1,num_regs] );
cell_ii     = nan( [1,num_regs] );

for ii = 1:num_regs


    [xx,yy] = getBBpad( props(ii).BoundingBox, imsize, 3 );
    [XX_ii,YY_ii] = meshgrid(xx,yy);
    
    ims_ii  = image_med(yy,xx);
    imd_ii  = imaged(yy,xx);
    
    mask_ii = (fr_label(yy,xx)==ii);
    
    [max_med_ii(ii), ind] = max(  imd_ii(mask_ii) );
    
    tmp = imd_ii(mask_ii);
    max_disk_ii(ii) = tmp(ind);
   
    tmp = XX_ii(mask_ii);
    x_ii(ii) = tmp(ind);
   
    tmp = YY_ii(mask_ii);
    y_ii(ii) = tmp(ind);
   
    x_ = x_ii(ii);
    y_ = y_ii(ii);
    %% assign to the correct cell
    
    %% figure out which cell
    ss_dist = [numel(yy),numel(xx)];

    cells_label = data.regs.regs_label(yy,xx);
    cells_mask = logical(cells_label);
    mask_tmp = zeros( ss_dist );
    
    mask_tmp( y_-yy(1)+1, x_-xx(1)+1 ) = 1;
    dist_tmp = bwdist( mask_tmp );

    list = cells_label(cells_mask);
    [~,ind_min] = min( dist_tmp(cells_mask) );
    cell_num = list(ind_min);

    %% Add the focus to the right cell
    if ~isempty( cell_num )
        cell_ii(ii) = cell_num;
        
        if debug_flag
           plot( x_, y_, '.r' );
           text( x_, y_, num2str( max_disk_ii(ii), '%1.2g' ));
        end
    end
    
    
end


%% Assign to cells.


for ii = 1:data.regs.num_regs

    list_ind = find(cell_ii == ii);
    
    max_med_ii_  = max_med_ii(list_ind);
    max_disk_ii_ = max_disk_ii(list_ind);
    x_ii_        = x_ii(list_ind);
    y_ii_        = y_ii(list_ind);
    cell_ii_     = cell_ii(list_ind);
    
    [max_disk_ii__,ord] = sort( max_disk_ii_, 'descend' );
    max_med_ii__        = max_med_ii_(ord);
    x_ii__              = x_ii_(ord);
    y_ii__              = y_ii_(ord);
    
    ind = find(max_disk_ii__ > 0.333*max_disk_ii__);
    if numel(ind) > CONST.trackLoci.numSpots(numc)
        ind = ind(1:CONST.trackLoci.numSpots(numc));
    end
    
    nfocus = numel(ind);

    
    max_disk_ii__ = max_disk_ii__(ind);
    max_med_ii__  = max_med_ii__(ind);
    x_ii__        = x_ii__(ind);
    y_ii__        = y_ii__(ind);
    
    %% make the subtracted mask;
    mask = data.CellA{ii}.mask;
    xx   = data.CellA{ii}.xx; 
    yy   = data.CellA{ii}.yy; 
    
    
    mask_mod = zeros( size( mask) );
    for jj = 1:nfocus
        
        yp = y_ii__(jj)+1-yy(1);
        xp = x_ii__(jj)+1-xx(1);

        if (xp>0) && (xp<numel(xx)+1) && (yp>0) && (yp<numel(yy)+1)
            mask_mod(yp, xp) = 1;
        end
    end
    mask_mod = imfilter( mask_mod, hd, 'replicate' );
    mask_mod = (mask_mod>0.5);
    
    mask_ii = and(mask,~mask_mod);
    im_ii = image_(yy,xx);
    
    Imean = mean( im_ii(mask_ii));
    Istd  = std(  im_ii(mask_ii));

    
    focus = focus0;

    for jj = 1:nfocus
        Iten = max_disk_ii__(jj)-A*Imean;
        score = Iten/(sqrt(A)*Istd);
        
        focus(jj).r               = [x_ii__(jj),y_ii__(jj)];
        focus(jj).score           = score;
        focus(jj).intensity_score = Iten;
        focus(jj).intensity       = Iten;
        focus(jj).b               = 1;
        focus(jj).error           = nan;
        
        %focus(jj).r = focus(jj).r + data.CellA{ii}.r_offset-[1,1];
        focus(jj).shortaxis = ...
            (focus(jj).r-data.CellA{ii}.coord.rcm)*data.CellA{ii}.coord.e2;
        focus(jj).longaxis = ...
            (focus(jj).r-data.CellA{ii}.coord.rcm)*data.CellA{ii}.coord.e1;
    end

    sc = [focus(:).score];
    focus = focus( ~isnan(sc) );
    
    data.CellA{ii} = setfield( data.CellA{ii}, fieldname, focus );

end


end
% 
% function junk
% %% Add the focus to the right cell
% if ~isempty( cell_num )
%     if isfield( data.CellA{cell_num}, fieldname )
%         focus_old = [getfield( data.CellA{cell_num}, fieldname ),...
%             focus];
%     else
%         focus_old = [focus];
%     end
%     
%     data.CellA{cell_num} = setfield( data.CellA{cell_num}, fieldname, focus_old );
% end
% 
% 
% %%
% 
% if disp_flag2
%     figure(1);
%     if cell_num
%    plot( x,y, '.b' ) 
%     else
%            plot( x,y, '.r' ) 
% 
%     end
%     
% end
% 
% end
% 
% function tmp
% %% Renormaze scores by cell
% % Use the std of the model that include subtracking foci
% for ii = 1:data.regs.num_regs
%     if ii == 207;
%        'hi'; 
%     end
%     if isfield( data.CellA{ii}, fieldname )
%         
%         xx = data.CellA{ii}.xx+crop;
%         yy = data.CellA{ii}.yy+crop;
%         
%         im_ii = model_(yy,xx);
%         
%         mask_ii = data.CellA{ii}.mask;
%         
%         dI_ii = std( double(im_ii(mask_ii)));
%         
%         focus = getfield( data.CellA{ii}, fieldname  );
%         num_loc = numel(focus);
%         
%         score_vec = [focus(:).score];
%         [~,ord_sort] = sort( score_vec, 'descend' );
%         focus = focus( ord_sort );
%         
%          for jj = 1:num_loc
%              focus(jj).intensity_score = focus(jj).intensity_score*dI0/dI_ii;
%              focus(jj).score           = focus(jj).score*dI0/dI_ii;
%              
%              focus(jj).shortaxis = ...
%                (focus(jj).r-data.CellA{ii}.coord.rcm)*data.CellA{ii}.coord.e2;
%              focus(jj).longaxis = ...
%                (focus(jj).r-data.CellA{ii}.coord.rcm)*data.CellA{ii}.coord.e1;
%          end
%         
%         
%     else
%         focus = focus0([]);
%     end
%     
%     data.CellA{ii} = setfield( data.CellA{ii}, fieldname, focus );
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function C = do_the_fit( ddd )
%         
%         gg = g_fit_fun(xx2,yy2,ddd(1),ddd(2),ddd(3),I0,ddd(4));
%         
%         %tmp = (im(iy2,ix2)-gg).*mask_crop(iy2,ix2);
%         %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
%         %tmp = (double(imm)-gg).*mask_crop(iy2,ix2);
%         
%         tmp = (double(imm)-gg);%.*mask_crop_for_fit(iy2,ix2);
%         C = sum(tmp(logical(mkk)).^2);
%         
%         if disp_flag
%             ddd(4);
%             
%             figure(3);
%             clf;
%             imshow( cat(3,ag([gg,imm]),ag([~mkk,~mkk]),0*[mkk,mkk]));
%             '';
%         end
%     end
% 
% end
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%
% %
% % Gaussian Fit function
% %
% %%%%%%%%%%%%%%%%%%
% 
% function C = g_fit_fun(X,Y,x0,y0,IG,I0,b1)
% 
% C = I0 + IG*exp( -((X-x0).^2+(Y-y0).^2)/(2*b1^2) );
% 
% end
% 




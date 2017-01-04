%% Model of cytoplasmic fluor in cell. Fit cell by cell.
% numc is channel number
% function adds the field cytoX, where cyto is the global cytofluor model 
% for channel X = numc.
function [ data ] = intCytoFluorFit( data, CONST, numc )

disp_flag = false;


persistent gf
persistent gf2
persistent opt_cyto

if isempty( gf )
    gf  = fspecial( 'gaussian', 21, 3 );
    gf2 = fspecial( 'gaussian', 21, 2 );
    opt_cyto =  optimset('MaxIter',25,'Display','off', 'TolX', 1/10);
end

opt_cyto.MaxIter = 100;






I0 = getfield( data, ['fl',num2str(numc),'bg'] );



    
    
image0 = double(getfield( data, ['fluor',num2str(numc)] ));

image = image0;
image(image<I0) = I0;
image = image-I0;

cyto  = zeros( size( image ) );
cyto_mult = zeros( [1,data.regs.num_regs] );

%% Loop through region by region
for ii = 1:data.regs.num_regs
% Now fit cell-by-cell

 % get the working region out here.
 [xx,yy] =  getBB(data.regs.props(ii).BoundingBox);

    image_ii = image( yy, xx );
    mask_ii  = (data.regs.regs_label(yy,xx)==ii);
    %cyto0_ii = cyto0(yy,xx );
    
    % this is the model distribution for cyto fluor in the cell
    cyto0_ii =  imfilter( sqrt(double(bwdist( ~mask_ii ))), gf2, 'replicate' );

    cyto_ii  = cyto(yy,xx );
    

    ddd = [max(image_ii(mask_ii))/max(cyto0_ii(mask_ii))];

    [ddd] = fminsearch( @do_the_fit_cyto,ddd,opt_cyto);
    
    
    cyto_ii = cyto_ii + ddd(1)*sqrt(double(bwdist( ~mask_ii )));
    
    if all(isfinite( cyto_ii(:)))
        cyto(yy,xx) = cyto_ii;
    end
end

% smooth the cyto image
cyto = imfilter( cyto, gf2, 'replicate' )+I0;

% cyto is the global cytofluor model for channel numc
data = setfield( data, ['cyto',num2str(numc)], cyto );


    %% fit function
    function C = do_the_fit_cyto( ddd )
        
        model_ii = intMakeBacker( cyto0_ii, [ddd,0]);
        
        dimage = model_ii-image_ii;
        
        C = sum(dimage(mask_ii).^2);
        
         if disp_flag
         
                  figure(3);
                  imshow( cat(3,ag([model_ii,image_ii]),ag([~mask_ii,~mask_ii]),0*[mask_ii,mask_ii]));
                  '';
                  
                  drawnow;
        end
        
    end
end


function im_backer = intMakeBacker( backer, ddd )

im_backer = abs(ddd(1)) * backer + ddd(2);

end

%% 
function cyto = intMakeCyto( regs_label, cyto_mult, gf2 )



[xx,yy] =  getBB(data.regs.props(ii).BoundingBox);



cyto0_ii =  imfilter( sqrt(double(bwdist( ~mask_ii ))), gf2, 'replicate' );


end






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


%% fit the focus positions, region by region.
for ii = 1:num_regs

    
    mask_g = and( ~data.regs(ii).fluor_label, data.regs(ii).mask );
    xx = data.regs(ii).xx;
    yy = data.regs(ii).yy;

    im = fluor(yy,xx);

    data.regs(ii).I0_meanS(born) = mean( double(im(mask_g)) );
    data.regs(ii).I0_stdS( born) = std(  double(im(mask_g)) );

    
    % get the fit for this frame, region ii
    
   %% Take care of kymo
 kymo_tmp = imrotate( fluor_proc(yy,xx), -data_mask.regs.props(ii).Orientation, 'bilinear'  );
 
 I0g = data.regs(ii).I0_meanS(born );
 kymo_tmp(kymo_tmp<I0g) = I0g;
 kymo_tmp = kymo_tmp - I0g;
 
 kymo_tmp = kymo_tmp.*double(data.regs(ii).kymo_mask);
 data.regs(ii).kymo(:,born) = sum( kymo_tmp );
    
end

end







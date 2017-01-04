function [im_det] = detFilter( im )

im = double(im);
x = -4:4;

[X,Y] = meshgrid( x, x);
R2 = X.^2+Y.^2;

b = 1.5;


gau = 1/(2*pi*b) * exp( -R2/(2*b) );

f_xx = ((X/b).^2-1/b).*gau;
f_yy = ((Y/b).^2-1/b).*gau;

f_xy = X.*Y.*gau/b^2;


im_xx = imfilter( im, f_xx, 'replicate' );
im_yy = imfilter( im, f_yy, 'replicate' );
im_xy = imfilter( im, f_xy, 'replicate' );

im_det = im_xx.*im_yy-im_xy.^2;
im_det(im_det<0) = 0;

im_det = im_det.^.5;

end


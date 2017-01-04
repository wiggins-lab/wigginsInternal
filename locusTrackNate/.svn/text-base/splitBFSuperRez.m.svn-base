function [] = copyBF( BF_image, image_number )

%% THIS FUNCTION IS USED TO SPLIT A DUAL-VIEW BRIGHTFIELD IMAGE FROM THE
% SUPER-SCOPE INTO TWO IMAGES, SAVE THEM IN A FORMAT THAT IS READABLE BY
% superseggeropti.m AND MOVE THEM INTO SEPERATE DIRECTORIES. 

mkdir 'bf1';
mkdir 'bf2';

im = imread(BF_image);

im1 = im(:,1:256) ;
im2 = im(:,257:512) ;

for ii = 1:9
    
    imwrite(im1,['bf1/split_data_t00',num2str(ii),'c1.tif'],'tif','Compression','none');
    imwrite(im2,['bf2/split_data_t00',num2str(ii),'c1.tif'],'tif','Compression','none');    
end

for ii = 10:99
    
    imwrite(im1,['bf1/split_data_t0',num2str(ii),'c1.tif'],'tif','Compression','none');
    imwrite(im2,['bf2/split_data_t0',num2str(ii),'c1.tif'],'tif','Compression','none');    
    
end

if image_number >=100

for ii = 100:image_number
    
    imwrite(im1,['bf1/split_data_t',num2str(ii),'c1.tif'],'tif','Compression','none');
    imwrite(im2,['bf2/split_data_t',num2str(ii),'c1.tif'],'tif','Compression','none');    
    
end

end

end
function [] = splitSuperRez(dirname)


%% THIS FUNCTION IS USED TO SPLIT A DUAL-VIEW FLUORECENT IMAGE FROM THE
% SUPER-SCOPE INTO GREEN AND RED IMAGES, AND SAVE THE OUTPUT TO FILENAMES
% THAT ARE INTERPRETABLE BY superseggeropti.m 

if ~exist([dirname,'split1']);
mkdir 'split1';
end

if ~exist([dirname,'split2']);
mkdir 'split2';
end

%%%% 2-8-2011 NOTE : IN THE CURRENT SETUP IM1 IS THE mCHERRY CHANNEL %%%%
%%%% AND IM2 IS THE GFP CHANNEL (I THINK)

file_filter  = '*.tif';

contents = dir (  file_filter );
num_im = numel( contents );


for ii = 1:9
     
    jj = ii-1;
    name = ['img_00000000',num2str(jj),'_Sequence_000.tif'] ;

    im = imread(name) ;
    
    im1 = im(:,1:256);
    im2 = im(:,257:512) ;
        
    imwrite(im1,['split1/', 'split_data_t00',num2str(ii),'c2.tif'],'tif','Compression','none') ;
    imwrite(im2,['split2/', 'split_data_t00',num2str(ii),'c3.tif'],'tif','Compression','none') ;
    
end

    name = ['img_00000000',num2str(9),'_Sequence_000.tif'] ;
    im = imread(name) ;
    
    im1 = im(:,1:256);
    im2 = im(:,257:512) ;
        
    imwrite(im1,['split1/', 'split_data_t0',num2str(10),'c2.tif'],'tif','Compression','none') ;
    imwrite(im2,['split2/', 'split_data_t0',num2str(10),'c3.tif'],'tif','Compression','none') ;

for ii = 11:99
     
    jj = ii-1;
    name = ['img_0000000',num2str(jj),'_Sequence_000.tif'] ;

    im = imread(name) ;
    
    im1 = im(:,1:256);
    im2 = im(:,257:512) ;
    
    imwrite(im1,['split1/', 'split_data_t0',num2str(ii),'c2.tif'],'tif','Compression','none') ;
    imwrite(im2,['split2/', 'split_data_t0',num2str(ii),'c3.tif'],'tif','Compression','none') ;
    
end

    name = ['img_0000000',num2str(99),'_Sequence_000.tif'] ;
    im = imread(name) ;
    
    im1 = im(:,1:256);
    im2 = im(:,257:512) ;
        
    imwrite(im1,['split1/', 'split_data_t',num2str(100),'c2.tif'],'tif','Compression','none') ;
    imwrite(im2,['split2/', 'split_data_t',num2str(100),'c3.tif'],'tif','Compression','none') ;


if (num_im >=101)
    
    for ii = 101:num_im
     
        jj = ii-1;
        name = ['img_000000',num2str(jj),'_Sequence_000.tif'] ;

        im = imread(name) ; 
    
        im1 = im(:,1:256);
        im2 = im(:,257:512) ;
    
    imwrite(im1,['split1/', 'split_data_t',num2str(ii),'c2.tif'],'tif','Compression','none') ;
    imwrite(im2,['split2/', 'split_data_t',num2str(ii),'c3.tif'],'tif','Compression','none') ;
    
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BatchJackie(dirname_,skip,clean_flag,res)
%
% This function processes a raw data set by
% (1) Aligning and cropping time series
% (2) Organizing files into directories
%        xy points are put in their own dir's
%        phase, fluor files put in their own dir's
%        makes seg and cell directories
% (3) Segmenting the frames into cell regions
% (4) Linking the regions between time steps
% (5) finding loci and calculating fluor statisticsdata
% (6) Putting complete cells into the cell dir.
%
% dirname_ : dir containing raw tif files
%
% skip     : The segmentation is performed every skip files
%          : this skip is very useful for high frequency timelapse where
%          : cells would switch back and forth between one and two
%          : segments leading to errors that are difficult to resolve.
%          : This segments every skip files, then copies the segments into
%          : the intermediate frames.
%
% clean_flag : Set this to be true to start from scratch and reseg all the
%            : files regardless of whether any seg files exist. If this
%            : flag is false, use existing segments, if they exist and
%            : new segs if they don't yet exist.
%
% res      : is a string that is passed to loadConstants(Mine).m to load
%          : the right constants for processing.
%
% Paul Wiggins, 4/17/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trackOptiSumLoci(dirname, CONST, header)

if ~exist('header')
    header = [];
end



dirseperator = filesep;
if(nargin<1 || isempty(dirname))
    
    dirname = '.';
    dirname = [dirname,dirseperator];
else
    if dirname(length(dirname))~=dirseperator
        dirname = [dirname,dirseperator];
    end
end


dirname_seg    = [dirname,'seg',dirseperator];


% Get the track file names...
contents=dir([dirname_seg '*_err.mat']);
if isempty( contents )
    clist.data = [];
    clist.def={};
    clist.gate=[];
else
    data_c = loaderInternal([dirname_seg,contents(1).name]);
    
end


% make num_t
dirname_phase  = [dirname,'phase',dirseperator];
contents = dir( [dirname_phase  '*.tif']);
num_t = numel( contents );


% reset values for nc
contents = dir([dirname,'fluor*']);
num_dir_tmp = numel(contents);
nc = 1;
num_c = 1;

for i = 1:num_dir_tmp
    if (contents(i).isdir) && (numel(contents(i).name) > numel('fluor'))
        num_c = num_c+1;
        nc = [nc, str2num(contents(i).name(numel('fluor')+1:end))+1];
    end
end



% 
% % crop box
% if exist( [dirname,filesep,'raw_im'] ,'dir')
%     disp('BSSO: aligned images exist');
%     if exist([dirname_,filesep,'raw_im',filesep,'cropbox.mat'],'file')
%         tmp = load( [dirname_,filesep,'raw_im',filesep,'cropbox.mat'] );
%         crop_box_array = tmp.crop_box_array;
%     else
%         crop_box_array = cell(1,10000);
%     end
% end
% 

        
        
%% generate sum image
%loads all fluorescence channels
contents_nc = cell([1,nc-1]);

for j = 2:num_c
    contents_nc{j-1} = dir( [dirname,'/','fluor',num2str(j-1),'/*.tif'] );
end


%%
sum_im = [];
for i = 1:num_t
    
    %i
    for j = 2:num_c
        
        name = contents_nc{j-1}(i).name;
        
        fluor_tmp(:,:,j-1) = imread( [dirname,'/','fluor',num2str(j-1),'/',name] );
        
    end
    
    if isempty( sum_im );
        sum_im = double( fluor_tmp );
    else  
        sum_im = sum_im +  double( fluor_tmp );
    end
    
end



%% generate fluor regions
disp_flag = true;

for mm = 1:num_c - 1
    data(mm).sum_im    = sum_im(:,:,mm);
    data(mm).mask_cell = data_c.mask_cell; 
    data(mm).mask_bg   = data_c.mask_bg; 
    data(mm).regs      = data_c.regs;
    data(mm).phase     = data_c.phase;
    data(mm).num_t     = num_t;
    
    data_mm = makeFluorRegions(data(mm), CONST,dirname, disp_flag );
    data(mm).im_smooth = data_mm.im_smooth;
    data(mm).mask_fluors = data_mm.mask_fluors;
    data(mm).fluor_label = data_mm.fluor_label;
    data(mm).fluor_props = data_mm.fluor_props;
end

clear data_mm;

%% Fit locii 
opt     =  optimset('MaxIter',100,'Display','off');
 

%setting up array which remembers which cells have had foci
%data(1) because regs.num_regs is number of cells, so same for all
%fluorescence channels - just pick one data
with_intensity = zeros(num_c - 1, data(1).regs.num_regs);

%records which frames are measured as having intensity
framesWithIntensity = zeros(length(with_intensity),num_t);

h = waitbar( 0, 'Wow! there are only...' );
for i = 1:num_t
    
    waitbar( i/num_t, h );
    
    
    
    
    for j = 2:num_c
        
        name = contents_nc{j-1}(i).name;

        
        fluor_tmp(:,:,j-1) = imread( [dirname,'/','fluor',num2str(j-1),'/',name] );
    end
    
    for mm = 1:num_c - 1
        disp_flag = false;
        %tic;
        [data_mm,with_intensity(mm,:),framesWithIntensity] = ...
            findFocusT6(framesWithIntensity, with_intensity(mm,:), fluor_tmp(:,:,mm), ...
            data(mm), CONST, opt, disp_flag, i );
        data(mm).cell = data_mm.cell;
        data(mm).framesWithIntensity = framesWithIntensity;
        clear data_mm;
    end
    %toc
    
    %saves data every 25 frames so can start to look at data
%     if mod(i,25) == 0
%         save( [dirname_xy,'/seg/data.mat'], 'data', '-v7.3' );
%     end
    
end
close(h);


% array of the numbers of cells which have foci
for mm = 1: num_c - 1
    with_intensity_mm = with_intensity(mm,:);
    with_intensity_mm(with_intensity_mm == 0) = [];
    data(mm).with_intensity = with_intensity_mm;
end

disp('Starting to write cell files');

with_intensity_unique = unique(with_intensity);
with_intensity_unique(with_intensity_unique == 0) = [];

recordNeighborTracesT6(dirname, data, with_intensity_unique, framesWithIntensity);

disp('cell files recorded');

save( [dirname,'/seg/data.mat'], 'data', '-v7.3' );
     

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function data = loaderInternal( filename );
data = load( filename );
end



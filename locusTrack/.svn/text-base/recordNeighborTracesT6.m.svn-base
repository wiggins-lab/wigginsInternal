%%this function writes cell files for cells which are specified by
%%with_intensity - in my case, with_intensity refers to cells which have shown
%%noticeable foci at any frame in this series

%in case of multiple fluorescence channels, if property is same for all
%channels (i.e. BoundingBox), by default use data from first fluorescence
%channel

function [] = recordNeighborTracesT6(dirname_xy, data, with_intensity, framesWithIntensity)

for jj = 1:length(with_intensity)
    
    kk = with_intensity(jj);
    
    %record which times have focus for each cell
    timesWithIntensity = [];
    for time = 1:size(framesWithIntensity,2)
        if framesWithIntensity(kk,time) == kk
            timesWithIntensity = [timesWithIntensity, time];
        end
    end
    currentcell.times_with_intensity = timesWithIntensity;
        
    
    %gives coordinates of bounding box around this cell
    ss = size(data(1).sum_im);
    
    [xx_cell,yy_cell] = getBBpad( data(1).regs.props(kk).BoundingBox, ss, 8 );
    
    regs_label_cell = data(1).regs.regs_label;
    
    regs_label_cell = unique(regs_label_cell(yy_cell(1):yy_cell(end),xx_cell(1):xx_cell(end)));
    regs_label_cell(regs_label_cell==0) = []; %remove zeros in regs_label
    regs_label_cell(regs_label_cell==kk) = []; %remove own identity from regs_label
    
    
    currentcell.regs_label = kk;
    
    for qq = 1:length(data)
        for pp = 1:length(data(qq).cell(kk).trace)
            currentcell.fluor(qq).trace(pp).I = data(qq).cell(kk).trace(pp).I;
        end
    end
    
    %add traces of neighbors
    for nn = 1:length(regs_label_cell)
        currentcell.neighbor(nn).regs_label = regs_label_cell(nn);
        
        %add all fluorescence traces of neighbors
        for qq = 1:length(data)
            for pp = 1:length(data(qq).cell(regs_label_cell(nn)).trace)
                currentcell.neighbor(nn).fluor(qq).trace(pp).I = data(qq).cell(regs_label_cell(nn)).trace(pp).I;
            end
        end
    end
    
    
    %load in cropped fluorescent images
    for qq = 1:length(data)
        fluorimages_qq = dir(fullfile([dirname_xy,'fluor',num2str(qq),'/*.tif']));
%     end    
%     for qq = 1:length(data)
        for ii = 1:length(fluorimages_qq) %ii is index of frame, not index of cell
            str = ['fluor',num2str(qq)];
            im = imread(fullfile(dirname_xy,str,fluorimages_qq(ii).name));
%             currentcell.fluor(qq).frame(ii) = im(yy_cell, xx_cell,:);
            cell_im_fluor( 1:length(yy_cell),1:length(xx_cell),ii) = im(yy_cell, xx_cell, :);
        end
        
        currentcell.fluor(qq).frame = cell_im_fluor;
    end
    
    %we're only working with one brightfield image, so only load one in
    
    phaseimages = dir(fullfile(dirname_xy,'phase/*.tif'));
    
    ii = length(phaseimages); %ii is index of frame, not index of cell
    im = imread(fullfile(dirname_xy,'phase',phaseimages(ii).name));
    cell_im_phase = im(yy_cell, xx_cell, :);
    
    
    currentcell.phase = cell_im_phase;
    
    newname = sprintf('cell%05d',kk);
    save([dirname_xy,'/cell/',newname,'.mat'],'-STRUCT','currentcell');
    
    %make movie for each cell
    
%     nFrame = numel(currentcell.fluor1);
%     
%     for ii = 1:nFrame
%         
%         back = data.phase{1};
%         epp = 0.33;
%         
%         imcolor = cat(3, epp*ag( back ), ...
%             ag(data.fluor1{ii})+epp*ag( back ), ...
%             epp*ag(back ));
%         imshow( imcolor );
%         
%         
%         M(ii) = getframe;
%         
%     end
%     
%     movie2avi(M, [dirname_xy,'/cell/',newname,'.avi']);
    
    clear currentcell;
    
end
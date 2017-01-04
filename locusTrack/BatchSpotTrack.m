function BatchSpotTrack(dirname_)
%batch function of running eSpot on several directories

if exist('loadConstantsMine');
    loadConstantsMine
else
    loadConstants
end



if ~exist('disp_flag') || isempty(disp_flag);
    disp_flag = 0;
end

if ~exist('file_filter') || isempty(file_filter);
    file_filter = '*.tif';
end

dirseperator = filesep;
% cmpstr=computer;%code for determining which computer this is running under
% if cmpstr(1:3)=='PCW'
%     dirseperator='\';%for pc
% elseif cmpstr(1:3)=='GLN'
%     dirseperator='/';%for linux
% elseif cmpstr(1:3)=='MAC'
%     dirseperator='/';%for mac (cheack dir seperator)
% end

if(nargin<1 || isempty(dirname_))
    dirname=['.',dirseperator];
else
    if dirname_(length(dirname_))~=dirseperator
        dirname=[dirname_,dirseperator,dirseperator];
    else
        dirname=dirname_;
    end
end


contents=dir([dirname 'cell*.mat' ]);



h = waitbar( 0, 'Track Spots');

num_cells = numel(contents);

for i=1:num_cells
    %i = 38
    
    
    dataname=[dirname,contents(i).name];
    data = load(dataname);
    
    
    num_im = numel(data.CellA);
    
    for j = 1:num_im
        
        waitbar((i-1+j/num_im)/num_cells,h,['Track Spots--Cell: ', ...
            num2str(i),'/',num2str(num_cells),' Frame: ',...
            num2str(j),'/',num2str(num_im)]);
        
        %j = 30
        celld = data.CellA{j};
        
        if isfield(celld, 'fluor1');
            celld.im_spot_1 = celld.fluor1;
        end
        
        if isfield(celld, 'fluor2');
            celld.im_spot_2 = celld.fluor2;
        end
        
        %celld.r_offset = celld.BB(1:2);
        celld.id       = data.ID;
        celld.A        = sum(double(logical(celld.mask(:))));
        celld.im_back  = celld.phase;
        celld.cell_num = data.ID;
        
        
        
        
        %data.CellA{j} = SSOMakeCell( celld, 1 );
        data.CellA{j} = intTrackSpotsDev( celld );
               
        
    end
    
    
    
    save(dataname,'-STRUCT','data');
    
end

close(h);
end
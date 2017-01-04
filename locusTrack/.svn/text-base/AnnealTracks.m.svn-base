function [data] = AnnealTrack( xx, map, fl );

E0 = 1e9;

ss = size(xx);

nt = ss(1);
nx = ss(2);

llink = cell(1,nx);
E     = intE( xx, map, fl, llink );

Nt = 1000;

for t = 1:Nt
    
    T = tempSchedule(t,Nt)
    
    
    llink0 = llink;
    E0     = E
    % perturb configuration
    llink = stepLink( llink, map );
    E     = intE( xx, map, fl, llink );
    
    
    DE = E0-E
    
    if rand > exp( DE/T )
        llink = llink0;
        E = E0;
        disp('reject')
    else
        
    end
    

    intDrawTrack( xx, map, fl, llink )
    drawnow;
    '';
    
end



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% intPlotTrack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intDrawTrack( xx, map, fl, llink )

nlink = numel(llink);


    clf;
    hold on;

map_list = unique(map);
map_list0 = [];

for i = 1:nlink
    
    ntrack = numel(llink{i});
    
    for j = 1:ntrack
        intPlotTrack( xx, map, fl, llink{i}(j), i )
    end
    
    map_list0 = [map_list0, llink{i}];
end

map_list0 = map_list(~ismember( map_list, map_list0 ));

ntrack = numel(map_list0);

for j = 1:ntrack
    intPlotTrack( xx, map, fl, map_list0(j), 0 )
end



end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% intPlotTrack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intPlotTrack( xx, map, fl, i, ii )
global pixelsize

cc = {'r','g','c','b','m','k'};

if ~ii
    ii = numel(cc);
end

xx_ = xx;
xx_(~fl) = NaN;


map_max = max(map(:));
ss = size(xx);

tt    = (1:ss(1))'*ones(1,ss(2));
tt_   = tt';
tt_ = tt_(:);
map_ = map';
map_ = map_(:);
xx_  = xx_';
xx_  = xx_(:);

ind = find(map_(:)==i);

plot( tt_(ind), xx_(ind)*pixelsize+0.05*ii ,cc{ii});


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% intPlotTrack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = intE( xx, map, fl, llink )
ss = size(map);
nt = ss(1);
nx = ss(2);

mapp = map';
xxp  = xx';
mapp = mapp(:);
xxp  = xxp(:);

E = 0;

map_list0 = [];
map_list = unique(map);

Et0 = 1;
Nt0 = 0;

El0 = 1;
Nl0 = 0;

Eol = 20;
Nol = 0;
Ex =  4;

Etot = 0;
DX = 0;



nlink = numel(llink);
for i = 1:nlink
    
    lnn = zeros(ss);
    
    llink{i} = sort( llink{i} );
    ntrack = numel(llink{i});
    for j =1:ntrack
        mm = llink{i}(j);
        lnn(mm==map) = 1;
                
        if j~=ntrack
            ind1 = max(find(mm==mapp));
            mm = llink{i}(j+1);
            ind2 = min(find(mm==mapp));
            DX = DX + abs(xxp(ind1)-xxp(ind2));
            
        end
    end
    
    E = E + Eol*sum(sum(lnn')==2) + Et0*sum(sum(lnn')==0);
    
end


nn = zeros(ss);
map_list0 = map_list(~ismember( map_list, map_list0 ));
n_map0 = numel(map_list0);
for i = 1:n_map0
    mm = map_list0(i);
    nn(mm==map) = 1;
end

E = E + El0*sum(nn(:)) + Ex * DX;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function llink = stepLink( llink, map );

map_max = max(map(:));
nl = numel( llink );

lnum = floor(nl*rand + 1);
mnum = floor(map_max*rand + 1);

ind = ( mnum==llink{lnum} );

if any(ind)
    llink{lnum} = llink{lnum}(~ind);
else
    llink{lnum} = [llink{lnum},mnum];
end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    T = tempSchedule(t,Nt)

t0 = Nt*.6;
dt = Nt/5;
E0 = 10;

T = E0*exp( (t0-t)/dt );

end
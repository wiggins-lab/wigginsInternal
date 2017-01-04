%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = bestLinkTmp( xx, xx_ref)
MAX_JUMP_X = 3;
JUMP       = MAX_JUMP_X^2;

n1 = numel( xx );
n2 = numel( xx_ref );

if ~isempty( xx_ref )
    
    
    xx     = reshape(xx,1,n1);
    xx_ref = reshape(xx_ref,n2,1);
    
    XX     = ones(n2,1)*xx;
    XX_ref = xx_ref*ones(1,n1);
    
    dist2 = (XX-XX_ref).^2;
    
    [dist2s, ords] = sort(dist2);
    
    run_flag = 1;
    
    while run_flag
        map =    ords(1,:);
        mapd=    dist2s(1,:);
        %dist2s
        %ords
        
        if numel(map) == numel(unique(map))
            break;
        end
        
        
        [dist2s, ords] = intInc( dist2s, ords );
        
        
    end
    
    
    map(mapd>JUMP) = 0;
else
    
    map = zeros(1,n1);
end
end


function [dist2s, ords] = intInc( dist2s, ords )

map =    ords(1,:);
umap = unique(map);

n_umap = numel(umap);

for j = 1:n_umap
    ind = find(map==umap(j));
    
    if numel(ind) > 1
        
        [tmp, ind2] = min(dist2s(1,ind));
        
        ind = ind(ind~=ind(ind2));
        
        for k = ind
            dist2s(:,k) = dist2s( [2:end,1], k );
            ords(:,k) = ords( [2:end,1], k );
            
        end
    end
end
end
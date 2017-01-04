%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [ind_min_dist, distance^2 ] = bestLink2( x1, x2 )
%
% This function matches the x positions in x1 with the positions in x2 by
% finding the minimum distance between them. The output is the an array of 
% inices and distances. match(1) is the index of x2 which is the smallest
% distance to x(1). DX2s(1) is this distance squared.
%
% Paul Wiggins, 8/12/2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [match,DX2s] = bestLink2( x1, x2 )


nl = numel(x1);

x1 = reshape(x1, 1, nl);
x2 = reshape(x2, nl, 1);

XX1 = ones(nl,1)*x1;
XX2 = x2*ones(1,nl);

DX2 = (XX1-XX2).^2;

[DX2s,ord] = sort(DX2);

% if iter is on, it makes the map one-to-one. This option goes both ways...
% introducing problems as well as fixing them. For now I am switching it
% off but I am leaving the code here since I think when other parts of the
% tracking code our changed, we might be able to turn this function on.


%iter

match = ord(1,:);
DX2s  = DX2s(1,:);
match(isnan(DX2s)) = 0;

    function iter
        
%         DX2s
%         ord
         ind = unique( ord(1,:) );
%                 
%         '';

        
        if numel(ind) ~= nl
            for j = ind
                jj = find(j==ord(1,:));
                
                [tmp,ord2] = sort(DX2s(1,jj));
                jj = jj(ord2);
                
                for j_ = jj(2:end)
                    ord(:,j_) = ord([2:end,1],j_);
                    DX2s(:,j_) = DX2s([2:end,1],j_);
                    
                end
            end
            
            iter;
            
        end
    end


end
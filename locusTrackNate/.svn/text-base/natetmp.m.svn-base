%nate tmp
% THIS FUNCTION BUILDS A MATRIX OF DISPLACEMENT VECTORS FROM THE
% OUTPUT OF THE getLocusTracks_Dev.m function
%
% 2-23-2011 PAW and NJK (mostly PAW)

function [dx, dx_op_or, dx_op, x_op] = natetmp( data );

% figure(20);
% clf;
% figure(20);
% hold on

xxp = data.xx'; 

ll = 0*data.xx;
dx = ll*NaN;
ss = size(ll);
nl = zeros(1,ss(1));
x_sort = ll*NaN;

for ii = 1: ss(1);
    tmp   = unique( data.ltrack(ii,:) );
    ntmp = numel(tmp);
    ll(ii,1:ntmp) = tmp;
    nl(ii) = ntmp;
end


for ii = 2:ss(1); % for the length of the possible track
    
    if nl(ii) ~= nl(ii-1); % if the number of unique tracks changes, set the disp = NaN
        dx(ii-1,:) = NaN;
    else
        
        
        
        for jj = 1:nl(ii) % for the number of unique tracks
            
            ind_eq = find( ll(ii,jj) == data.ltrack(ii,:)); % find the pointer in ltrack
            ind_eq = ind_eq(1); % choose the row value (ostensibly the track it associates with)
            
            if (ll(ii,jj) ~= 0 && xxp(ll(ii,jj)) ~= xxp(data.ltrack(ii-1,ind_eq)))
                
                dx(ii,jj) = xxp(ll(ii,jj)) - xxp(data.ltrack(ii-1,ind_eq)); % calculate the displacement
                
                dx_op(ii,jj) = sign(xxp(ll(ii,jj))); % record the sign of the trajectory coordinate
                
                x_op(ii,jj) = xxp(ll(ii,jj));
                
            end
            
            
            %            xxp(ll(ii,jj))
            %            xxp(data.ltrack(ii-1,ind_eq))
            %
            %            data.xx(max(1,(ii-2)):(ii+0),:)
            %            data.ltrack(max(1,(ii-2)):(ii+0),:)
            %
            %
            %            'hi'
            %            if nl(ii) > 1
            %                 'hi'
            %            end
        end
        
        
        
        
        
    end
    
    
    
end


max_op = max(find(dx_op(:,1) == -dx_op(:,2))); % find the order of the tracks (which is positive and which is negative

dx_op_or = dx_op(max_op,:); % dx_op_or is a 2-col vector, either [-1 1] or
% [1 -1], giving the relative orientation
% of the split loci, such that the red and
% green tracks can be matched to each other


% PLOT FOR DEBUGGING
% num = sum(double(~isnan(dx))')';
% ind = find( num == 2);
%
% figure(20);
% plot(dx(ind,1),dx(ind,2),'o')
% axis square
% hold on

end

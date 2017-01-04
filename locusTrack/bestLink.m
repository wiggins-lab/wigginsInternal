%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% findCombo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match0 = bestLink( x1, x2 )

n1 = numel(x1);
x1 = reshape(x1,1,n1);

n2 = numel(x2);
x2 = reshape(x2,n2,1);

x1 = ones(n2,1)*x1;
x2 = x2*ones(1,n1);

dist2 = (x2-x1).^2;

    
nn = numel(dist2(1,:) )+1;

nnn1 = sum(~isnan(dist2(1,:)));
nnn2 = sum(~isnan(dist2(:,1)));

%nnn1 = numel(dist2(1,:));
%nnn2 = numel(dist2(:,1));


NNN1 = nnn1+1;
NNN2 = nnn2+1;

MAX_JUMP = 3;

JUMP = MAX_JUMP^2;

match = zeros(1,nn);

DDmin = [];

moveN( 1 )

    function moveN( n_where )
        %n_where
        %match
        if n_where == NNN1+1
            
            DD = 0;
            for i = 1:NNN1
                j = match(i);
                if j == 0;
                    dist2_ = JUMP;
                else
                    try
                        if i == NNN1
                            dist2_ = JUMP;
                        elseif j == NNN2
                            dist2_ = JUMP;
                        else
                            dist2_ = dist2(j,i);
                        end
                    catch
                        '';
                    end
                    if isnan(dist2_)
                        dist2_ = JUMP;
                    end
                end
                DD = DD + dist2_;
                
            end
            
            %             match
            %             DD
            if isempty(DDmin) || (DD < DDmin)
                DDmin = DD;
                match0 = match;
            end
            
        else
            
            
            for n_current = 0:NNN2
                
                if n_current == 0
                    match(n_where) = n_current;
                    moveN(n_where+1);
                elseif all( (n_current ~= match(1:n_where-1)) )
                    match(n_where) = n_current;
                    moveN(n_where+1);
                end
            end
        end
    end




end
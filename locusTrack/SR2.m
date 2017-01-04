function [ind1_no_match, ind1_match, indA_match] = SR2( x1, y1, xA, yA, r )





n1 = numel( x1 );
nA = numel( xA );

if (n1 > 0) && (nA > 0)

    try 
        XA = ones(n1,1)*xA;
        X1 = x1'*ones(1,nA);
    catch
       'hi' 
    end
        YA = ones(n1,1)*yA;
        Y1 = y1'*ones(1,nA);
        
        R = sqrt(((XA-X1).^2 + (YA-Y1).^2)/(r)^2);
        
        ind1_no_match    = find(~any(R < 1,2))';
        ind1_match       = find(any(R < 1,2))';
        [tmp,indA_match] = min(R(ind1_match,:)');

        ind_singles      = find(sum(R < 1,1)==1);        
        A = ismember( indA_match, ind_singles);

        indA_match = indA_match(A);
        ind1_match = ind1_match(A);
%        indA_nol = find(sum(double(R < 1),1)==1);       
%        ind1_match = indA_nol(ind1_match_);
        
        
       
end

if n1 == 0
    ind1_no_match = [];
    ind1_match    = [];
    indA_match    = [];
elseif nA == 0
    ind1_no_match = 1:n1;
    ind1_match    = [];
    indA_match    = [];
end

end

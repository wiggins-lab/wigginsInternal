ss = size(tra);
num_corr = 1;

for ii = 1:(ss(1)-num_corr)
    
    jj = ii+1;
    
    tra1(ii) = tra(jj,1) - tra(ii,1);
    tra2(ii) = tra(jj,2) - tra(ii,2);
    
end

figure(20)
plot(tra1,tra2)
    

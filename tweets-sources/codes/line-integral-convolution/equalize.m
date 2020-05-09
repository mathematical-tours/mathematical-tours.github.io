function F = equalize(F,u)

n = size(F,1);
[tmp,I] = sort(F(:));  
F(I) = sort(u(:)); % sort(randn(n*n,1)); % linspace(0,1,n*n);
    
end

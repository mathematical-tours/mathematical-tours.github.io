function [h,k] = Convergents(a)

% Convergents - the fraction is h./k

n = length(a);

h2 = 0; k2 = 1; 
h1 = 1; k1 = 0;
h = []; k = [];
for i=1:n
    [h(i),k(i)] = deal( a(i)*h1+h2,a(i)*k1+k2 );    
    [h1,k1,h2,k2] = deal(h(i), k(i), h1, k1);
end

end
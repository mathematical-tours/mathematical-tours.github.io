function c = coef_spline_1(ud)
K = length(ud); alpha = sqrt(3)-2; 
b1 = alpha; b0=-6*b1/(1-b1^2);
c1 = zeros(K,1); c2 = zeros(K,1);
for i=2:K
    c1(i) = ud(i)+b1*c1(i-1);    
end
for i=(K-1):-1:1
    c2(i) = ud(i)+b1*c2(i+1);    
end
c = b0*(c1+c2-ud);
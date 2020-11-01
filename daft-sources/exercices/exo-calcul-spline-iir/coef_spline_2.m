function c = coef_spline_2(ud)
K = length(ud); alpha = sqrt(3)-2; 
c = zeros(K,1); d = zeros(K,1);
for i=2:K
    d(i) = 6*ud(i)-alpha*d(i-1);    
end
c(K) = -6*alpha/(1-alpha^2)*(2*d(K)-6*ud(K));
for i=(K-1):-1:1
    c(i) = alpha*( c(i+1)-d(i+1) );    
end
function Bc = brownian_bridge(P,K,a,b,sigma)

t = linspace(0,1,P)';
% generate a bunch of bridges
B = (randn(P-1,K)+1i*randn(P-1,K))/sqrt(P);
B = [zeros(1,K); cumsum(B)];
B = B - t * B(end,:);
B = B*sigma;
Bc = B + repmat(1-t,[1 K])*a + repmat(t,[1 K])*b;

end
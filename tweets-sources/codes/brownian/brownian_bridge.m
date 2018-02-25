function Bc = brownian_bridge(P,K,sigma,a,b)

t = linspace(0,1,P)';
% generate a bunch of bridges
B = (randn(P-1,K)+1i*randn(P-1,K))/sqrt(P);
B = [zeros(1,K); cumsum(B)];
if nargin>3 && not(isempty(a))
    B = B - t * B(end,:);
    B = B*sigma;
    Bc = B + repmat(1-t,[1 K])*a + repmat(t,[1 K])*b;
else
    Bc = B*sigma;
end

end
% Transformée de Walsh rapide
function y = fwt(x)
N = length(x);  % N doit être une puissance de 2
if(N==1) y = x; return; end;
P = N/2;
x = [fwt(x(1:P)) ; fwt(x((P+1):N))];
y = zeros(N,1);
y(1:P) = x(1:P) + x((P+1):N);
y((P+1):N) = x(1:P) - x((P+1):N);
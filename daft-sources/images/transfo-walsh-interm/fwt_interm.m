% Transformée de Walsh intermédiaire rapide
function y = fwt_interm(x, alpha)
N = length(x);  % N doit être une puissance de 2
if(N==1) y = x; return; end;
P = N/2;
x = [fwt_interm(x(1:P), alpha) ; fwt_interm(x((P+1):N), alpha)];
y = zeros(N,1);
y(1:P) =     sqrt(2)*( cos(alpha)*x(1:P) + sin(alpha)*x((P+1):N) );
y((P+1):N) = sqrt(2)*( sin(alpha)*x(1:P) - cos(alpha)*x((P+1):N) );
function y = czt_simpson(x,alpha)
N = length(x); N0 = (N-1)/2; y = zeros(N,1);
y(1) = x(1)/3; y(N) = x(N)/3;
y( 2*(1:(N0-1))+1 ) = 2/3*x( 2*(1:(N0-1))+1  );
y( 2*(0:(N0-1))+2 ) = 4/3*x( 2*(0:(N0-1))+2  );
y = czt(y,alpha);
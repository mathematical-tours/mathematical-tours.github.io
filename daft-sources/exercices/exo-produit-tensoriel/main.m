N = 128;

x = zeros(N,1);
x(100) = 1;
y = trigo_tensoriel(x, pi/10);

    plot(0:N-1,y,'k.:');
    axis square;
%    axis([0,N-1,-100,100]);
    axis tight;

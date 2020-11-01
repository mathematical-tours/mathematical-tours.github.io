global u1;
global u2;
global N;

N = 32;

f = zeros(N,1);
f(1) = 1;
f(2) = 1;
f(3) = 1;
%f(4) = 1;

x = (0:N-1)*2*pi/N;
% f = sin(x)';

u1 = sqrt(N)*f+fht(f);
u2 = sqrt(N)*f-fht(f);

x = (-N/2+1):N/2;
c = [((N/2+1):N)';(1:N/2)'];

% subplot(1,2, 1);
% plot(x, f(c), 'k*-' );
% legend('f');
% axis([-N/2+1,N/2,-0.05,1.05]);
% title('Signal d''origine');

subplot(1,2, 1);
plot( x, u1(c), 'k.:' );
title('U_+(f)');
axis tight;
axis square;

subplot(1,2, 2);
plot( x, u2(c), 'k.:' );
title('U_-(f)');
axis tight;
axis square;

% saveas(gcf, '../racine-carree-tfd', 'eps');
% saveas(gcf, '../racine-carree-tfd', 'png');

% pause;
clf;
nb = 5;
t = 0;
for l = 0:1/(nb-1):1
    t = t+1;
    subplot(1,nb, t);
    y = interm(l);
    y = y/max(abs(y));
    plot( x, real( y(c) ), 'k:.' ); % , x, imag( interm(l) ), 'k+-' );
    text( -12,0.9, ['\', sprintf( 'lambda=%0.1f', l )] ); 
    axis([-N/2+1,N/2,-0.3,1]);
    axis square;
end


saveas(gcf, '../transfo-hartley-interm', 'eps');
saveas(gcf, '../transfo-hartley-interm', 'png');

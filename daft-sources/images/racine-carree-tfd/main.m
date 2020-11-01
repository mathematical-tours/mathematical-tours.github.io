global u1;
global u2;
global u3;
global u4;
global N;

N = 20;

f = zeros(N,1);
f(1) = 1;
f(2) = 1;
f(3) = 1;
%f(4) = 1;

rev = [1;(N:-1:2)'];
fs = ( f+f(rev) )/2;
fa = ( f-f(rev) )/2;


u1 = real( sqrt(N)*fs+fft(fs) );
u2 = real( sqrt(N)*fs-fft(fs) );
u3 = real( sqrt(N)*fa+i*fft(fa) );
u4 = real( sqrt(N)*fa-i*fft(fa) );

x = (-N/2+1):N/2;
c = [((N/2+1):N)';(1:N/2)'];

% subplot(1,2, 1);
% plot(x, f(c), 'k*-' );
% legend('f');
% axis([-N/2+1,N/2,-0.05,1.05]);
% title('Signal d''origine');

subplot(1,2, 1);
plot( x, u1(c), 'k*:', x, u2(c), 'k+:' );
legend('U_+(f)', 'U_-(f)');
axis tight;
axis square;
title('Vecteurs propres pour \pm 1');

subplot(1,2, 2);
plot( x, u3(c), 'k*:', x, u4(c), 'k+:' );
legend('V_+(f)', 'V_-(f)');
axis tight;
axis square;
title('Vecteurs propres pour \pm i');

saveas(gcf, '../racine-carree-tfd', 'eps');
saveas(gcf, '../racine-carree-tfd', 'png');
pause;
clf;
nb = 5;
t = 0;
for l = 0:1/(nb-1):1
    t = t+1;
    subplot(1,nb, t);
    y = interm(l);
    y = y/max(abs(y));
    plot( x, real( y(c) ), 'k:.' ); % , x, imag( interm(l) ), 'k+-' );
    text( -8.5,0.9, ['\', sprintf( 'lambda=%0.1f', l )] ); 
    axis([-N/2+1,N/2,-0.05,1]);
    axis square;
end


saveas(gcf, '../racine-interm-tfd', 'eps');
saveas(gcf, '../racine-interm-tfd', 'png');

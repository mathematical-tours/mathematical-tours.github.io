% convolution et multiplication de polynomes

N = 11;
P = 6;
% 1+2X+X^3-X^4+X^5
f = zeros(N,1);
f(1)=1;
f(2)=2;
f(3)=0;
f(4)=1;
f(5)=-1;
f(6)=1;


% X-X^2+2X^3+2X^5
g = zeros(N,1);
g(1)=0;
g(2)=1;
g(3)=-1;
g(4)=2;
g(5)=0;
g(6)=2;

% X+X^2+5X^4+8X^6+4X^8-3X^7-2X^9+2X^10
c = 1:N; % [(N/2+1:N):N,1:N/2]';
h = real( ifft(fft(f(c)).*fft(g(c))) );
h = h(c);

subplot(1,3,1);
plot(0:N-1, f, 'k*');
axis tight;
axis square;
title('1+2X+X^3-X^4+X^5');


subplot(1,3,2);
plot(0:N-1, g, 'k*');
axis tight;
axis square;
title('X-X^2+2X^3+2X^5');


subplot(1,3,3);
plot(0:N-1, h, 'k*');
axis tight;
axis square;
title('X+X^2+5X^4+8X^6+4X^8-3X^7-2X^9+2X^{10}');

saveas(gcf, '../convolution-polynome', 'eps')
saveas(gcf, '../convolution-polynome', 'png')
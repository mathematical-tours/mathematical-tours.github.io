	clear;
clf;

N = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noyeau de la gaussienne.
et = 0.15;
% calcule d'une gaussienne recentrée
f = zeros(N,N);
for i=(-N/2+1):N/2
for j=(-N/2+1):N/2
	x = pi*i/N;
	y = pi*j/N;
	f(i+N/2,j+N/2) = exp( -(x^2+y^2)/et );
end
end
s = sum(sum(f));
f = f/s;	% renormalisation


% LA fonction de matlab.
P = 40;
g = peaks(P);
g = g + 0.8*rand(P,P);

subplot(1,3,1);
xx = (-P/2+1):P/2;
surf(xx, xx, g);
axis([-P/2+1 P/2  -P/2+1 P/2 -5,6]);
axis square;
title('f');

subplot(1,3,2);
xx = (-N/2+1):N/2;
surf(xx, xx, f);
axis([-N/2+1,N/2, -N/2+1,N/2, 0,1/s]);
axis square;
title('g');

% convole les deux fonction.
h = conv2(g,f, 'same');
subplot(1,3,3);
surf(h);
axis([1, P, 1, P, -5,6]);
axis square;
title('f*g');


saveas(gcf, '../convolution-2d', 'eps');
saveas(gcf, '../convolution-2d', 'png');
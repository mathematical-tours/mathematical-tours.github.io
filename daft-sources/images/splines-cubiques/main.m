%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N+2 points
N = 100;
h = 1/(N+1);
xx = 0:h:1;
yy = rand(N+2,1);

% les dérivées extrèmes
der_0 = 0;
der_1 = 0;
N = length(xx);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saisie des points
% axes('position', [0,0,1,1]); 
% [xx,yy] = ginput;
yy = zeros(8,1);
yy(1) = 0.1;
yy(2) = 1;
yy(3) = 0.5;
yy(4) = 1.3;
yy(5) = 0.2;
yy(6) = 0.8;
yy(7) = 0.8;
yy(8) = 1.3;
xx = (1:8)';

h = 1/(N-1);
xx_t = min(xx):h:(max(xx)+h);


%% spline libre
yy_t = SPLINE(xx', [der_0 yy' der_1] ,xx_t);
subplot(1,2,1);
plot(xx_t-1, yy_t, xx-1, yy, 'o');
axis([-0.5,7.5,0,1.5]);
axis square;
title('Spline libre');


%% spline nonlibre
yy_t = SPLINE(xx', yy' ,xx_t);
subplot(1,2,2);
plot(xx_t-1, yy_t, xx-1, yy, 'o');
axis([-0.5,7.5,0,1.5]);
axis square;
title('Spline "not-a-knot"');


saveas(gcf, '../spline-cubique', 'eps');
saveas(gcf, '../spline-cubique', 'png');

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20, 'XTick', [], 'YTick', []);


% #samples
n = 3000; p = n;

name = 'chirp-sine'; 
name = 'warped-sine';

switch name
    case 'chirp-sine'
        x = chirp(0:n-1,0,n,1/100);
        y = cos(2*pi*5*(0:p-1)/p);
    case 'warped-sine'
        f = 12; % frequency
        psi = @(t).5 + .5*sign(2*t-1) .* abs(2*t-1).^3;
        x = cos(2*pi*f*psi( (0:n-1)/n ));
        y = cos(2*pi*f*(0:p-1)/p);
        plot(x);
end

[d,IX,IY] = dtw(x,y, 'euclidean');


% Display distance matrix
D = abs( repmat(x(:), [1 p]) - repmat(y(:)', [n 1]) );

clf; hold on;
imagesc((1:p)/p,(1:n)/n,D); 
plot(IY/p,IX/n, 'g', 'LineWidth', 3); 
axis image; axis off;
colormap gray(256);
saveas(gcf, [rep name 'match.png']);

h = 1.5;
tx = linspace(0,1,n); 
ty = linspace(0,1,p);

clf; plot(tx,x, 'r', 'LineWidth', 2);
SetAR(1/2);
saveas(gcf, [rep name '-1.eps'], 'epsc');
clf; plot(ty,y, 'b', 'LineWidth', 2);
SetAR(1/2);
saveas(gcf, [rep name '-2.eps'], 'epsc');

return;

clf; hold on;
plot(tx,x-h, 'b');
plot(ty,y+h, 'r');
sk = 20; 
for s=1:sk:length(IX)
    i = IX(s); j = IY(s);
    %plot([tx(i) ty(j)],[x(i)-h, y(j)+h], 'k');
    plot([tx(i) ty(j)],[-h, +h], 'k');
end
% plot(linspace(0,1,p),y+h, 'r');
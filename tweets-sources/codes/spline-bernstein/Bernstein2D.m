%%
% In 2D

addpath('../toolbox/');
rep = MkResRep();

m = 1024;
x = linspace(0,1,m)';

f = @(x)(.2+x)/1.2.*exp(2*2i*pi*x);

%
q=70;
nlist = unique(1+round(linspace(0,1,q).^7 * 300));
%
clf; hold on;
plot(f(x), 'k', 'LineWidth', 2);
for it=1:length(nlist)
    n = nlist(it);
    s = (it-1)/(length(nlist)-1);
    y = BernsteinBasis(x,n+1) * f( (0:n)'/n );
    %
    plot(y, 'r', 'Color', [s 0 1-s]);
    axis equal; axis([-1 1 -1 1]); 
    box on; set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    drawnow;
end

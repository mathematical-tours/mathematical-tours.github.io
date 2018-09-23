%%
% Display of evolution of Alpha Shapes


addpath('../toolbox/');
rep = MkResRep();

if not(exist('x'))
n = 200; 
x = rand(n,1);
y = rand(n,1);
%
x = []; y = [];
clf; hold on;
while true
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    x(end+1) = a; 
    y(end+1) = b;
end
x = x(:); y = y(:);

end
n = length(x);

shp = alphaShape(x,y);
alist = alphaSpectrum(shp); 
alist = alist(1:2:end);
alist = alist(end:-1:1);

q = 100;
alist = linspace(0,.2,q);

%q = 100;
%alist = linspace(0,.2,q);


alist = alphaSpectrum(shp); 
alist = alist(end:-1:1);
alist = alist(alist<.2);

q = length(alist);
for i=1:q
    t = (i-1)/(q-1);
    col = [t 0 1-t];
    shp.Alpha = alist(i);
    clf; hold on;
    h = plot(shp);
    set(h,'FaceColor',.1*col+.9);
    set(h,'EdgeColor',col);
    set(h,'LineWidth',2);
    plot(x,y, 'k.', 'MarkerSize', 20);
    axis equal; axis([0 1 0 1]); axis off;
    drawnow;
    saveas(gcf, [rep 'evol-' znum2str(i,3) '.png']);
end

% AutoCrop(rep, 'evol-');
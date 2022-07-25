addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


% input (d+1) x d matrix
simplex = randn(3,2);

f0 = @(x,y)sqrt(x.^2+y.^2);

% Himmelblau's function
f0 = @(x,y)(x.^2+y-11).^2 + (x+y.^2-7).^2;
f0 = @(x,y)f0(x*6,y*6);

f = @(x)f0(x(1),x(2));
clf;
[x_opt, n_feval, Xlist]  = nelder_mead( simplex, f, 1 );

n = 256;
t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
R = f0(X,Y);
R = rescale(R);

[~,I] = sort(R(:));
R(I) = linspace(0,1,n*n);


clf; hold on;
imagesc(t,t,R');
contour(t,t,R',linspace(0,1,r), 'k');
colormap(parula(r-1));
caxis([0 1]);
axis image; axis off;
it = 0;
for i=1:length(Xlist)
    s = (i-1)/(length(Xlist)-1); s = s^.5;
    x = Xlist{i};
    plot(x([1:3 1],1), x([1:3 1],2), '.-', 'MarkerSize', 20, 'LineWidth', 2, 'color', [s 0 1-s]);
    % h = area(x([1:3 1],1), x([1:3 1],2), 'FaceColor', 'r', 'LineStyle', 'None');
    % h.FaceAlpha = 0.5;
    axis([-1 1 -1 1]);
    drawnow;
    for j=1:3
        it = it+1;
        mysaveas(it);
    end
end




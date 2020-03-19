n = 512; 


addpath('../toolbox/');
rep = MkResRep();

t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
% contour( X.*Y./(X.^2+Y.^2) );

F = X.*Y./(X.^2+Y.^2);
F = (X.^2) .* Y .* sqrt(X.^2+Y.^2) ./ (X.^4+Y.^2);
F = (X.^2) .* Y ./ (X.^2+Y.^2);

q = 50;
for it=1:q
    s = (it-1)/q;
    clf;
    surf(F);
    axis tight;
    shading interp;
    % camlight;
    camlight;
    view(s*180,40);
    axis off
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png'], 'png');
end
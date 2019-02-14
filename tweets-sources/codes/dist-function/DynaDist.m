%%
% dynamic distfunc.

addpath('../toolbox/');
rep = MkResRep('dynamic');

subdivide = @(f,h)cconvol( upsampling(f), h);
h = [1 4 6 4 1]; % cubic B-spline
h = 2*h/sum(h(:));

% grid
n = 256;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

%%
% First we create a 2D closes polygon.

[f0,f1,f2] = selectPoly();


q = 50; 

for it=1:q
    s = (it-1)/(q-1);
    if s<1/3
        s1 = 3*s;
        f = (1-s1)*f0 + s1*f1;
    elseif s<2/3
        s1 = 3*(s-1/3);
        f = (1-s1)*f1 + s1*f2;
    else
        s1 = 3*(s-2/3);
        f = (1-s1)*f2 + s1*f0; 
    end
    % subdivide
    g = f;
    for k=1:6
        g = subdivide(g,h);
    end
    % distance
    D = distmat( [real(g) imag(g)]',[X(:)';Y(:)'] );
    D = reshape(min(D,[],1),[n n]);
    %%%% draw %%%%
    % draw
    r = 15;
    clf; hold on;
    imagesc(t,t,D');
    colormap(parula(r+1));
    %
    plot(real(g([1:end 1])), imag(g([1:end 1])), 'r', 'LineWidth', 2);
    %
    axis equal;
    axis([0 1 0 1]);
    axis off;
    drawnow;
    saveas(gcf, [rep 'evol-' znum2str(it,2) '.png'], 'png');
end

% AutoCrop(rep, 'evol');



% solve for tsp from an animation.


addpath('./toolbox-lsap/');
addpath('../toolbox/');
rep = MkResRep();

n = 1000;
xy = rand(n,2);

% generate data through error diffusion

name = 'turing';
name = 'lovelace';

names = {'turing' 'lovelace'};

method = 'simple';
method = 'floyd';


xy = {};
for k=1:length(names)
    n = 128*4;
    f = load_image(names{k},n);
    f = sum(f,3);
    [~,I] = sort(f(:)); f(I) = linspace(0,1,n*n);
    
    % target points
    h = 1-f;
    p = 2*5000;
    g = ErrorDiffusion(h/sum(h(:)) * p );
    %
    s = round(.1*n):round(.9*n);
    g = g(s,s);
    % g(:,1) = 0; g(:,end) = 0; g(1,:) = 0; g(end,:) = 0;
    %
    [y,x] = ind2sub(size(g), find(g(:)));
    x = (x-1)/(size(g,1)-1); y = (y-1)/(size(g,1)-1);
    y = 1-y;
    
    clf;
    plot(x,y, '.');
    axis xy;    
    xy{k} = [x,y];    
end

p = min(size(xy{1},1), size(xy{2},1));
for k=1:2
    I = randperm(size(xy{k},1));
    I = I(1:(size(xy{k},1)-p));
    xy{k}(I,:) = [];
end

% comptute OT.

a = ones(p,1) / p;
C = distmat(xy{1}',xy{2}').^2;
C1 = int32( round(C*1e6) );
[J,varrho,u,v] = hungarianLSAP(C1);



q = 50; % #frames
tlist = linspace(0,1,q);
for k=1:length(tlist)
    t=tlist(k);
    Xt = (1-t)*xy{1} + t*xy{2}(J,:);
    %
    clf; hold on;
    plot(Xt(I,1), Xt(I,2), '.', 'color', [1-t;0;t]);
    axis equal; axis([0 1 0 1]);  box on;
    set(gca, 'XTick', [], 'YTick', []); drawnow;
    drawnow;
    saveas(gcf, [rep 'ot-' znum2str(k,2) '.png']);
end


niter = 20000;
[~,I] = solveTSP( xy{1}, false, niter );
niter = 1000;


q = 50; % #frames
tlist = linspace(0,1,q);
for k=1:length(tlist)
    t=tlist(k);
    Xt = (1-t)*xy{1} + t*xy{2}(J,:);
    [~,I] = solveTSP( Xt, false, niter, I );
    %
    clf; hold on;
    plot(Xt(I,1), Xt(I,2), 'color', [1-t;0;t], 'LineWidth', 1);
    axis equal; axis([0 1 0 1]);  box on;
    set(gca, 'XTick', [], 'YTick', []); drawnow;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(k,2) '.png']);
end


% AutoCrop(rep, ['anim-'])

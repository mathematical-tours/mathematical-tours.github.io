%%
% Display evolution of marginal/conditional where the coupling is computed
% using Sinkhorn.

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
SetTickOff = @()set(gca, 'XTick',[], 'YTick',[]);
SetTickOn = @()set(gca, 'XTick',[0 1/2 1], 'YTick',[0 1/2 1]);

N = 256*2;
t = linspace(0,1,N)';
gauss = @(m,s)exp(-(t-m).^2/(2*s).^2);
normalize = @(x)x/sum(x(:));

s = .6;
itr = @(u,v,t)s * ( (1-t/.43)*u+t/.43*v );
A = @(t)gauss(.3*s+t,.05*s) + .5*gauss(.6*s+t,itr(.2,.05,t));
B = @(t).5*gauss(1-s*.2-t,.04*s) + .8*gauss(1-s*.5-t,itr(.04,.15,t)) + .5*gauss(1-s*.8-t,.04*s);
vmin = .01;
A = @(t)normalize(vmin + A(t));
B = @(t)normalize(vmin + B(t));

if 0
for it=1:q
    s = (it-1)/(q-1);
    a = A(s*.43);
    b = B(s*.43);
    clf; hold on;
    h = area(t, a, 'FaceColor', 'r', 'EdgeColor', 'r');
    h.FaceAlpha = 0.5;
    h = area(t, b, 'FaceColor', 'b', 'EdgeColor', 'b');
    h.FaceAlpha = 0.5;
    axis tight; SetAR(1/2); SetTickOff();     
    drawnow;
end
end


epsilon = .015; % for sinkhorn
% square cost function
C = abs(t-t').^2;

q = 50;

for it=1:q

    s = (it-1)/(q-1);
    a = A(s*.43);
    b = B(s*.43);
    
    options.niter = 500;
    options.tol = 1e-12;
    [P,f,g,Err] = sinkhorn(C,a,b,epsilon,options);
    
    
    clf;
    area(t, a, 'FaceColor', 'r', 'EdgeColor', 'r');
    axis tight; SetAR(1/2); SetTickOff();
    saveas(gca, [rep 'a-' znum2str(it,2) '.png']);
    
    clf;
    area(t, b, 'FaceColor', 'b', 'EdgeColor', 'b');
    axis tight; SetAR(1/2); SetTickOff();
    saveas(gca, [rep 'b-' znum2str(it,2) '.png']);
    
    r = 16; % #levellines
    L = linspace(0,1,r); L(end)=[];
    m = linspace(0,1,r-1)';
    
    S = P; S = S/max(S(:));
    clf; hold on;
    imagesc(t,t,1-S');
    contour(t,t,1-S',L, 'k');
    colormap( (1-m)*[.5 0 .5] + m*[1 1 1] );
    caxis([0 1]);
    axis image; box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gca, [rep 'joint-' znum2str(it,2) '.png']);
    
    % Conditional
    S = P ./ b'; S = S/max(S(:));
    clf; hold on;
    imagesc(t,t,1-S');
    contour(t,t,1-S',L, 'k');
    colormap( (1-m)*[1 0 0] + m*[1 1 1] );
    caxis([0 1]);
    axis image; box on;
    set(gca, 'XTick', [], 'YTick', []);
    saveas(gca, [rep 'cond-a-' znum2str(it,2) '.png']);
    
    
    % Conditional
    S = P ./ a; S = S/max(S(:));
    clf; hold on;
    imagesc(t,t,1-S');
    contour(t,t,1-S',L, 'k');
    colormap( (1-m)*[0 0 1] + m*[1 1 1] );
    caxis([0 1]);
    axis image; box on;
    set(gca, 'XTick', [], 'YTick', []);
    saveas(gca, [rep 'cond-b-' znum2str(it,2) '.png']);


end

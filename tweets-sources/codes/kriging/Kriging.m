%%
% Display kriging, aka RBF with confidence intervals.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

p = 8; 
n = 1024; % for vizualization
t = linspace(0,1,n)';

rand('state', 123134);
x1 = rescale(cumsum(1+rand(p,1)*10), .05, .95);
x2 = rescale(cumsum(1+rand(p,1)*10), .05, .95);
f1 = rescale(rand(p,1), .1,.9);
f2 = rescale(rand(p,1), .1,.9);

sc = .05; % scale
K = @(x)exp(-x.^2/(2*sc^2));

q = 70; 

for it=1:q
    s = (it-1)/(q-1);
    x = (1-s)*x1 + s*x2;
    f = (1-s)*f1 + s*f2;
    % mean
    w = K(x-x')\f;
    F = K(t-x') * w;    
    % cov
    C = 1 - sum( ( inv( K(x-x') ) * K(x-t') ) .*  K(x-t') );
    S = .2*sqrt(max(C,0))';
    %
    lineProps.col = {[s 0 1-s]};
    clf; hold on; 
    mseb(t',F',S', lineProps);
    plot(x, f, '.', 'color', [s 0 1-s], 'MarkerSize', 25);
    % plot(t, F, 'color', [s 0 1-s], 'LineWidth', 2);
    axis([0 1 0 1]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end
  
%%
% Display of Brownian motion.

addpath('../toolbox/');
rep = MkResRep();

P = 1024*2; % #points
K = 10;  % #trajectory
sigma = 1;
alph = .1;
lw = 2;

%%
% brownian motion

Bc = brownian_bridge(P,1,sigma);
clf;
plot_colored(Bc,alph,lw);
r = 1.5;
axis equal; axis([-r r -r r]);
axis off;
saveas(gcf, [rep 'brownian.png'], 'png');

%%
% Brownian bridge

slist = [0 .05 .1 .2 .5];
for i=1:length(slist)
    sigma = slist(i);
    a = 0; b = 1;
    Bc = brownian_bridge(P,K,sigma,a,b);
    clf;
    plot_colored(Bc,alph,lw);
    r = 1.5;
    axis equal; axis([-.1 1.1 -.4 .4]);
    axis off;
    saveas(gcf, [rep 'bridge-' num2str(round(100*sigma)) '.png'], 'png');
end
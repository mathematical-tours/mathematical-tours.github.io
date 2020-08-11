%%
% Display the evolution of a Brownian motion.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);



n = 2048*4; % #sample for the simulation
k = 10; % number of paths

sigma = .5;
W = sigma*(randn(n,1)+1i*randn(n,1))/sqrt(n);
W = W-mean(W)+1/n;
B = [0;cumsum(W)];

q = 100;
plist = round( linspace(10,400,q) );
plist = round(2.^(2:.25:13)');

clf; hold on;
for it=1:length(plist)
    s = (it-1)/(length(plist)-1);
    p = plist(it);    
    I = ceil( p*(1/n:1/n:1) );
    W1 = W;
    for i=1:p
        W1(I==i) = mean(W(I==i));
    end
    B1 = [0;cumsum(W1)];
%    plot(B, 'b');
    clf;
    plot(B1, '-', 'color', [s 0 1-s], 'LineWidth', 2);
    axis equal;
    axis([-.1 1.1 -.6 .6]); axis off;
    drawnow;
    mysaveas(it)
end
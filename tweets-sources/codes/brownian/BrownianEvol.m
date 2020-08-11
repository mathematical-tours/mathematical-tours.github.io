%%
% Display the evolution of a Brownian motion.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 2048*8; % #sample for the simulation
k = 50; % number of paths

W = (randn(n,k)+1i*randn(n,k))/sqrt(n);
B = [zeros(1,k);cumsum(W)];

q = 200; 

col = distinguishable_colors(k);
xmax = 2;
for it=1:q
    s = (it-1)/(q-1);
    i = round( s*n + 1 );
    clf;  hold on;
    for m=1:k
        plot(B(1:i,m), 'color', .3*col(m,:)+.7, 'LineWidth', 1);
    end
    for m=1:k
        plot(real(B(i,m)),imag(B(i,m)), '.', 'color', col(m,:), 'MarkerSize', 20);
    end
    axis equal; 
    axis([-1 1 -1 1]*xmax); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('particles', it);
end

% display the evolution of the density
x = linspace(-xmax,xmax,512)';
for it=1:q
    s = min(1e-5 + (it-1)/(q-1),1); % variance
    h = exp(-(x.^2+x'.^2)/(2*s));
    if 0
        h1 = h/sum(h(:));
        sum(sum( h1 .* (x.^2) )) % should be ==s
    end
    r = 15; % #levelines
    clf; hold on;
    imagesc(x,x,h');
    contour(x,x,h',linspace(0,1,r), 'k');
    m = linspace(0,1,r-1)';
    colormap(m*[s 0 1-s] + (1-m)*[1 1 1]);
    caxis([0 1]);
    axis image; axis off;
    drawnow;
    mysaveas('density', it);
end
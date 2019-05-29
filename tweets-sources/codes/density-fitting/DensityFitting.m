%%
% Synthetic example for density fitting.

n = 40; % sample points.
N = 5*1e6; % parzen estimate

% remap to [0,1]
c = [0,3.5]; s = 6;
R = @(z).5+.5i + .5*(real(z)-c(1))/s + .5*1i*(imag(z)-c(2))/s;
R = @(z)imag(R(z)) + 1i*real(R(z));


f = @(z)2*sign(real(z)).*abs(real(z)).^1 + 1i * ( (imag(z)/9) - .6*real(z).^2 + 8.5 );

z = randn(n,1) + 1i*randn(n,1);
Z = randn(N,1) + 1i*randn(N,1);

clf; hold on;
plot(R(z), 'r.', 'MarkerSize', 20);
plot(R(f(z)), 'b.', 'MarkerSize', 20);
axis equal; axis([0,1,0,1]); box on; set(gca, 'XTick', [], 'YTick', []);

q = 40;
for it=1:q
    t = (it-1)/(q-1);
    zt = R( (1-t)*z + t*f(z) );
    Zt = R( (1-t)*Z + t*f(Z) );

    s = 6; % parzen smoothing 
    p = 400; % #pixels
    h = parzen2d(real(Zt),imag(Zt),p,s);
    h = h/max(h(:));
    u = linspace(0,1,p);


    r = 10; % #levellines
    m = linspace(0,1,r-1)';
    CM = m*[1-t 0 t] + (1-m)*[1 1 1];

    
    clf; hold on;
    % display quantized colormap
    lvl = linspace(0,1,r); lvl = lvl(2:end);
    clf; hold on;
    imagesc(u,u,h');
    contour(u,u,h',lvl, 'k');
    colormap(CM);
    caxis([0 1]);
    %
    % plot(Zt, 'b.', 'MarkerSize', 20);
    % plot(R(f(z)), 'r.', 'MarkerSize', 20);
    axis equal; axis([0,1,0,1]); box on; set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end


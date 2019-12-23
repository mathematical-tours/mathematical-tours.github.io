%%
% Show the equivalence sum / convolution

addpath('../toolbox/');
rep = MkResRep();


n = 256; % grid size
p = 2*1000; % #samples

t = [0:n/2,-n/2+1:-1]' / n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));



x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);
r = .2;
U = abs(X-1/2)<r & abs(Y-1/2)<r;

X = .5 + 1i*.5 + (2*rand(p,1)-1)*r + 1i * (2*rand(p,1)-1)*r;
Z = randn(p,1) + 1i * randn(p,1);

% convolution
r = 15; % #levellines
q = 50;
slist = linspace(1e-5,.1,q);

for it=1:q
    sigma = slist(it);
    s = (it-1)/(q-1);
    R = GFilt(U,sigma); R = R/max(R(:));
    %
    m = linspace(0,1,r-1)';
    CM = m*[s 0 1-s] + (1-m)*[1 1 1];
    %  
    clf; hold on;
    imagesc(x,x,R');
    A = linspace(0,1,r); A(1) = [];
    contour(x,x,R',A, 'k');
    colormap(CM);
    caxis([0 1]);
    axis image; axis off;
    saveas(gcf, [rep 'density-' znum2str(it,2) '.png']);
    
    % point cloud
    clf;
    plot(X + sigma * Z, '.', 'MarkerSize', 7, 'Color', [s 0 1-s]);
    axis equal; axis([0 1 0 1]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'points-' znum2str(it,2) '.png']);

end
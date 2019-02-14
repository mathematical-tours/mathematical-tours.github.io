addpath('../toolbox/');
rep = MkResRep();

name = 'hibiscus';

n = 512/2; 

f = load_image(name,n);
f = clamp(f/255);
d = size(f,3);

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

% compute the stack
p = 30*2; % #convolution images
smax = 50/2;
slist = linspace(1e-5,smax,p);
F = [];
for i=1:p
    for k=1:d
        F(:,:,k,i) = GFilt(f(:,:,k),slist(i));
    end
end

% compute image mask
t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
Gauss = @(m,s)exp(-((X-m(1)).^2+(Y-m(2)).^2)/(2*s^2));

% animation
q = 50; 
% animate
s = .3*2;

for it=1:q
    t = (it-1)/q;
    m = .6 * [cos(2*pi*t), sin(2*pi*t)];
    %
    M = round((1-Gauss(m,s)) * (p-1))+1;
    M = repmat(M,[1 1 d]);
    % selector
    [I,J,K] = meshgrid(1:n,1:n,1:k);
    g = F( I + (J-1)*n + (K-1)*n*n + (M-1)*n*n*d );
    imwrite(g, [rep 'foveated-' znum2str(it,2) '.png']);
    %
    clf
    imageplot(g);
    drawnow;
    % display level
    if 0
    R = Gauss(m,s);
    r = 15; % #levellines
    clf; hold on;
    imagesc(R);
    contour(R,linspace(0,1,r), 'k');
    colormap(parula(r-1));
    caxis([0 1]);
    axis image; axis ij; axis off;
    saveas(gcf, [rep 'map-' znum2str(it,2) '.png']);
    end
end

% AutoCrop(rep, 'map');


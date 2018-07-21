%%
% Test for Wasserstein flow of MMD losses in 1D.

addpath('toolbox/');
addpath('img/');

KernelType = 'variance'; 
KernelType = 'doublewell'; 

name = 'cross';
name = '2disk';
name = 'annulus';

rep = ['results/2d/' name '/' KernelType '/']; 
[~,~] = mkdir(rep);

%%
% The kernel


% dimension of the points
d = 2;

% x is of size (n,d)
Delta = @(x,z) repmat( reshape(x,[size(x,1) 1 d]), [1 size(z,1) 1]) ...
             - repmat( reshape(z,[1 size(z,1) d]), [size(x,1) 1 1]);
         
switch KernelType
    case 'doublewell'
        r0 = .1;
        phi = @(r)(r.^2-r0^2).^2 / 4; phi1 = @(r)r.*(r.^2-r0^2);
    case 'variance'
        phi = @(r)r.^2/2; phi1 = @(r)r;
end

% Gaussian kernel
Kn = @(x,z)sqrt(sum(Delta(x,z).^2,3))+1e-6;
K = @(x,z)phi( Kn(x,z) );
% its derivative with respect to first variable
Kt = @(x,z)phi1( Kn(x,z) );
K1 = @(x,z)Delta(x,z) .* repmat( Kt(x,z) ./ Kn(x,z), [1 1 d] );

% step size
tau = 30; niter = 100;



% load input measures
n = 150;
n0 = 60; % pixel in images
[Y,X] = meshgrid(linspace(0,1,n0), linspace(0,1,n0));
f = rescale(sum(load_image(name, n0),3), 0,1);
I = find(f<.5); % I = I(randperm(length(I))); I = I(1:n); I = I(:);
R = [X(I) Y(I)];
a = ones(length(I),1)/length(I);


% MMD distance and its gradient
dotp = @(x,y)sum(x(:).*y());
Energy  = @(X) dotp(K(X,X)*a,a)/2 - dotp(K(X,Y)*b,a) + dotp(K(Y,Y)*b,b)/2;
Grad = @(X) repmat(a,[1 d]) .* tensor_mult(K1(X,X),a);


% display the inputs
clf;
plot(R(:,1), R(:,2), '.', 'MarkerSize', 20);
axis equal; axis([0 1 0 1]); axis off;



normalize = @(x)x/sum(x);
% for the kernel estimator
s = .03;
q = 300; % image size for rendering

r = 15; % #levellines
t = linspace(0,1,q);

hx = KernelDensity2D(R,q,s);
clf; hold on;
imagesc(t,t,hx');
contour(t,t,hx',linspace(-1e-2,max(hx(:)),r), 'k');
colormap(parula(r-1));
caxis([0 max(hx(:))]);
axis image; axis off;
drawnow;
saveas(gcf, [rep 'input.png']);


Q = 50; % #display
ndisp = round(linspace(1,niter,Q));
k = 1;
       
X = R; 
for i=1:niter
    u = (i-1)/(niter-1);
    % display
    if i==ndisp(k)
        hx = KernelDensity2D(X,q,s);
        clf; hold on;
        imagesc(t,t,hx');
        contour(t,t,hx',linspace(-1e-2,max(hx(:)),r), 'k');
        colormap(parula(r-1));
        caxis([0 max(hx(:))]);
        axis image; axis off;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(k,2) '.png']);
        k = k+1;
    end
    % gradient step
    X = X - tau*Grad(X);
end


% AutoCrop(rep,'anim-')

%%
% Test for Wasserstein flow of MMD losses in 2D.

addpath('toolbox/');
addpath('img/');

name = 'variance'; 
name = 'aggregative';
name = 'repulsive'; 
name = 'attractive'; 


addpath('../toolbox/');
rep = MkResRep(name);

%%
% The kernel


% dimension of the points
d = 2;

% x is of size (n,d)
Delta = @(x,z) repmat( reshape(x,[size(x,1) 1 d]), [1 size(z,1) 1]) ...
             - repmat( reshape(z,[1 size(z,1) d]), [size(x,1) 1 1]);
Periodize = @(X)X.*(abs(X)<=1/2) + (X-1).*(X>1/2) + (X+1).*(X<-1/2);
Delta = @(x,z)Periodize(Delta(x,z));

% step size
tau = 30; niter = 400;

switch name
    case 'doublewell'
        r0 = .1;
        phi = @(r)(r.^2-r0^2).^2 / 4; phi1 = @(r)r.*(r.^2-r0^2);
    case 'variance'
        phi = @(r)r.^2/2; phi1 = @(r)r;
    case 'repulsive'
        m = .05;
        phi = @(r)1./(m+r); 
        phi1 = @(r)-1./(m+r).^2;
        tau = .2;
    case 'attractive'
        phi = @(r)r.^3/3; 
        phi1 = @(r)r.^2;
        tau = 10;
    case 'aggregative'
        m = .05;
        %
        psi  = @(r)r.^2 ./ (1+r).^3; 
        psi1 = @(r)2*r ./ (1+r).^3 - 3*r.^2 ./ (1+r).^4; 
        %
        phi = @(r)psi(r/m); 
        phi1 = @(r)psi1(r/m)/m;
        tau = 15*3;
end


%% 
% Display interaction potentials.

r = linspace(0,.4,1024);
clf;
plot(r, phi(r), 'k', 'LineWidth', 2); axis tight;
set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []); box on;
saveas(gcf, [rep 'potential.eps'], 'epsc');

% radial kernel
Kn = @(x,z)sqrt(sum(Delta(x,z).^2,3))+1e-6;
K = @(x,z)phi( Kn(x,z) );
% its derivative with respect to first variable
Kt = @(x,z)phi1( Kn(x,z) );
K1 = @(x,z)Delta(x,z) .* repmat( Kt(x,z) ./ Kn(x,z), [1 1 d] );

% load input measures
n = 600;
gauss = @(n,m,s)[randn(n,1)*s+m(1),randn(n,1)*s+m(2)]; 
R = [gauss(n/2,[.3 .3],.1);gauss(n/2,[.6 .7],.08)]; 

% uniform weights
n = size(R,1);
a = ones(n,1)/n;


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


if 0
hx = KernelDensity2D(R,q,s);
clf; hold on;
imagesc(t,t,hx');
contour(t,t,hx',linspace(-1e-2,max(hx(:)),r), 'k');
colormap(parula(r-1));
caxis([0 max(hx(:))]);
axis image; axis off;
drawnow;
saveas(gcf, [rep 'input.png']);
end

Q = 80; % #display
ndisp = round(linspace(1,niter,Q));
k = 1;

disp_type = 'density';
disp_type = 'points';
       
X = R; 
for i=1:niter
    u = (i-1)/(niter-1);
    % display
    if i==ndisp(k)
        clf; hold on;
        switch disp_type
            case 'density'
                hx = KernelDensity2D(X,q,s);
                imagesc(t,t,hx');
                contour(t,t,hx',linspace(-1e-2,max(hx(:)),r), 'k');
                colormap(parula(r-1));
                caxis([0 max(hx(:))]);
                axis image; axis off;
                drawnow;
            case 'points'
                plot(X(:,1), X(:,2), '.', 'MarkerSize', 25, 'Color', [u 0 1-u]);
                axis equal; axis([0 1 0 1]);
                set(gca, 'XTick', [], 'YTick', []); box on;
        end
        saveas(gcf, [rep 'anim-' znum2str(k,2) '.png']);
        k = k+1;
    end
	drawnow;
    % gradient step
    X = X - tau*Grad(X);
    X = mod(X,1);
end


% AutoCrop(rep,'anim-')

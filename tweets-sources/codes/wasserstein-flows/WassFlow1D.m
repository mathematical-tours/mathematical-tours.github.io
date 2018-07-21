%%
% Test for Wasserstein flow of MMD losses in 1D.

addpath('toolbox/');
 
EnergyType = 'variance';
EnergyType = 'doublewell'; 

rep = ['results/' EnergyType '/']; 
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20, 'XTick', [], 'YTick', []);


%%
% The kernel

Delta = @(x,z) repmat(x(:),[1 length(z)]) - repmat(z(:)',[length(x) 1]);
switch EnergyType
    case 'doublewell'
        r0 = .1;
        phi = @(r)(r.^2-r0^2).^2 / 4; phi1 = @(r)r.*(r.^2-r0^2);
    case 'variance'
        phi = @(r)r.^2/2; phi1 = @(r)r;
end

% step size
tau = 40; niter = 300;
% Gaussian kernel
K = @(x,z)phi(Delta(x,z));
% its derivative with respect to first variable
K1 = @(x,z)phi1(Delta(x,z));

% input measures
n = 2000; m = n;
a = ones(n,1)/n;
b = ones(m,1)/m;

% initial positions
X0 = [linspace(.2,.35,n/2), linspace(.7,.9,n/2)]';


% MMD distance and its gradient
dotp = @(x,y)sum(x(:).*y());
Energy  = @(X) dotp(K(X,X)*a,a)/2;
Grad = @(X) a .* ( K1(X,X)*a );

% for the kernel estimator
normalize = @(x)x/sum(x);
q = 500;
t = linspace(0,1,q);
mbar = @(hx,c)area(t, hx, 'EdgeColor', c, 'FaceColor', c);
s = .01;

hx0 = KernelDensity1D(X0,q,s);

vmax = max(hx0)*3;

Q = 50; % #display
ndisp = round(linspace(1,niter,Q));
k = 1;

X = X0; 
for i=1:niter
    u = (i-1)/(niter-1);
    % display
    if i==ndisp(k)
        hx = KernelDensity1D(X,q,s);
        clf; hold on;
        mbar(hx0,  .5 + .5*[1 0 0]);
        mbar(hx,  [1-u 0 u]);
        SetAR(1/2); box on;
        axis([0 1 0 vmax]);
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(k,2) '.png']);
        k = k+1;
    end
    % gradient step
    X = X - tau*Grad(X);
end


% AutoCrop(rep,'anim-')

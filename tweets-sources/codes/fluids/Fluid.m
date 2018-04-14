%%
% Fluid dynamics


name = 'random';
name = 'dirac';
name = 'random1';
name = 'dirac1';

addpath('../toolbox/');
rep = MkResRep(['dynamics/' name]);

n = 256/2;
normalize = @(V)V ./ repmat( max(1e-9,sqrt(sum(V.^2, 3))) , [1 1 2]);

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
G = @(s)exp(-t.^2/(2*s^2)); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

% grad and div 
a = [2:n 1]; b = [n 1:n-1]; % periodic bc
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f);
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) );
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) );
DeltaV = @(V)cat(3, Delta(V(:,:,1)), Delta(V(:,:,2)) );

% Fourier transform of Delta
[Y,X] = meshgrid(0:n-1,0:n-1);
mu = sin(X*pi()/n).^2; mu = -4*( mu+mu' );
% pseudo-inverse of Laplacian
mu(1) = 1; % avoid 0 division
DeltaInv = @(u)real( ifft2( fft2(u) ./ mu ) );
% Projection on incompressible flows.
ProjI = @(u)u + Grad(DeltaInv(Div(u)));



% Helper function: enforce periodicity.
periodic = @(P)cat(3, mod(P(:,:,1)-1,n)+1, mod(P(:,:,2)-1,n)+1 );
% Helper function: extend an image by 1 pixel to avoid boundary problems.
extend1 = @(f)[f f(:,1)];
extend = @(f)extend1(extend1(f)')';
% Helper function: bilinear interpolation on a grid.
myinterp = @(P1,f1,Pi)interp2( P1(:,:,2), P1(:,:,1),f1, Pi(:,:,2), Pi(:,:,1) );
% First we compute the initial and warped grids.
[Y,X] = meshgrid(1:n,1:n);  P = cat(3, X,Y);
[Y1,X1] = meshgrid(1:n+1,1:n+1); P1 = cat(3, X1,Y1);
% Defines the warping operator \(\Ww_U\).
W = @(f,U)myinterp( P1, extend(f), periodic( P - U ) );


% input image
imname = 'hibiscus';
f = rescale( sum(load_image(imname, n),3) );
% input field
s = 10;
V = cat(3, GFilt(randn(n),s), GFilt(randn(n),s) );
U = normalize(ProjI(V));

% example of advection along the flow
rho = .25;
niter = 12*4;
k = 0;
clf;
f1 = f;
for i=1:niter
    f1 = W(f1, rho*U );
    if mod(i,niter/4)==0
        k = k+1;
        imageplot(f1, strcat(['t=' num2str(i*rho)]), 2,2,k);
    end
end


% Set the viscosity \(\nu\) for the velocity field.
nu = 1/20;
% Use a larger viscosity \(\mu\) for the evolution of the density of particules.
mu = 1*nu;
% Extend the warping operator \(\Ww_U\) to work with vector fields as input.
Wt = @(V,U)cat(3, W(V(:,:,1),U), W(V(:,:,2),U) );
% We discretize the PDE's using some time step \(\tau\).
tau = .1;

% to display VF
sf = 1.2; % stretch factor
lw = 2;
q = 20;
k = round(n/q); % subsampling
plotvf = @(v,col)quiver(v(1:k:end,1:k:end,2),v(1:k:end,1:k:end,1), sf, 'Color', col, 'LineWidth', lw);

% initialization
switch name
    case 'random'
        s = 10;
        V = cat(3, GFilt(randn(n),s), GFilt(randn(n),s) );
        V = ProjI(V);
        g = f;
        %
        tau = .04;
        nu = 1/100;
        mu = 0;
        %
        niter = 300;
    case 'random1'
        s = 10;
        V = cat(3, GFilt(randn(n),s), GFilt(randn(n),s) );
        V = ProjI(V);
        %
        g = zeros(n);
        I = randperm(n*n); I = I(1:10); g(I) = 1;
        %
        tau = .04;
        nu = 1/100;
        mu = 0;
        %
        niter = 400;
    case 'dirac'
        V = zeros(n,n,2)+1e-3; 
        kappa = 50; 
        V(end/2,end/2,1)=-kappa;
        V(end/2,end/2,2)=kappa/2;
        V = ProjI(V);
        g = f;
        %
        tau = 5;
        nu = 1/100;
        mu = 0;
        %
        niter = 2000;
    case 'dirac1'
        V = zeros(n,n,2)+1e-3; 
        g = zeros(n);
        r = 30;
        g(end/2-r:end/2+r,end/2-r:end/2+r)=1;
        % g = double(randn(n)>0);
        kappa = 50; 
        % 
        % p = round([.7 .5]*n);
        % V(p(1),p(2),1)=-kappa; V(p(1),p(2),2)=kappa/2;
        %
        p = round([.1 .9]*n);
        % p = round([.5 .5]*n);
        V(p(1),p(2),1)=.5*kappa; V(p(1),p(2),2)=-.8*kappa;
        %
        V = ProjI(V);
        %
        % g = zeros(n);
        % I = randperm(n*n); I = I(1:100); g(I) = 1;
        %
        tau = 6;
        nu = 1/60;
        mu = 0;
        %
        niter = 3000;
end



ndisptot = 100;
ndisp = round(niter/ndisptot);


savemode = 1;
% 
clf; k = 0;
for i=1:niter
    % Advect
    g = W (g,tau*V);
    V = Wt(V,tau*V);
    % Diffuse
    V = V + tau*nu*DeltaV(V);
    g = g + tau*mu*Delta(g);
    % Project
    V = ProjI(V);
    % Display
    if mod(i,ndisp)==1
        k = k+1;
        t = (i-1)/(niter-1);
        % display field
        if savemode
        clf; hold on; 
        plotvf(V, [t 0 1-t]);
        axis equal; axis off;
        saveas(gcf, [rep 'flow-' znum2str(k,3) '.png'], 'png');
        end
        % display image
        clf; 
        imageplot(g);
        drawnow;
        if savemode
        imwrite(rescale(g), [rep  'img-' znum2str(k,3) '.png'], 'png');
        end
    end
end


% create gifs
% AutoCrop(rep, ['flow-'])
% Use % convert interp-*.png interp.gif to generate the gif using imagemagik




%%
% Folker-Planck on a grid.

addpath('../toolbox/');
rep = MkResRep('eulerian');

n = 256; % image size

% Laplacian
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per
% grad
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f)/2;
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) )/2;
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;



t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);


% driving force
U = cat(3,X,Y);

f0 = zeros(n,n);
c1 = [-.3 -.5]; c2 = [.5 .5];
f0 = (max(abs(X-c1(1)),abs(Y-c1(2)))<.4) + ...
     (max(abs(X-c2(1)),abs(Y-c2(2)))<.4);


 % heat, df/dt = Delta(f)
niter = 7000;
tau = .6;
kappa = 15;

q = 100;
k = 1;
ndisp = round(linspace(1,niter,q));
 
f  = f0;
for i=1:niter    
    if i==ndisp(k)
        u = (i-1)/(niter-1);
        % display 
        g = f/max(abs(f(:)));
        %
        r = 12; % #levellines
        clf; hold on;
        imagesc(t,t,g);
        contour(t,t,g,linspace(0,1,r), 'k');
        M = linspace(0,1,r-1)';
        colormap(M*[u 0 1-u] + (1-M)*[1 1 1]);
        caxis([0 1]);
        axis image; axis off; axis xy;
        drawnow;
        saveas(gcf, [rep 'eulerian-' znum2str(k,3) '.png'], 'png');
        k = k+1;
    end
    f = f - tau * Div( repmat(f, [1 1 2]) .* U / kappa + Grad(f) );
end


% AutoCrop(rep, 'eulerian-');
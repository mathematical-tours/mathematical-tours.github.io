%%
% Display of Laplacian spectrum on 2D domains.

name = 'france';
name = 'cat';
name = 'elephant';

addpath('../toolbox/');
rep = MkResRep(name);

n = 200;
f = load_image(name, n);
f = sum(f,3); f = double((f/max(f(:)))>.5);
if f(1)==0
    f = 1-f;
end
rho = f + 1e-5;

% grad operator with neumann BC
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);
dx = dx(2:end,:);
% with periodic BC
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);
dx(1,end) = -1;
%
Grad = [kron(dx,speye(n)); kron(speye(n),dx)];
%
Delta = Grad' * spdiags( [rho(:); rho(:)], 0, 2*n*n,2*n*n ) * Grad;

% do SVD
k = 40; % #eigenvectors
[U,S,V] = eigs(Delta,k+1, 'SM'); S = diag(S);
[~,I] = sort(S); U = U(:,I); U = U(:,2:end);

% display
r = 15; % #levellines
t = linspace(0,1,n);
for i=1:k
    u = reshape(U(:,i), [n,n]);
    u = u/max(abs(u(:)));
    u(f==1) = NaN;
    % display
    clf; hold on;
    imagesc(t,t,u);
    contour(t,t,u,linspace(-1,1,r), 'k');
    colormap(jet(r-1));
    caxis([-1 1]);
    axis image; axis off; axis ij;
    drawnow;
    pause(0.1);
    saveas(gcf, [rep 'eig-' znum2str(i,2) '.png' ]);
end

 % AutoCrop(rep, 'eig');

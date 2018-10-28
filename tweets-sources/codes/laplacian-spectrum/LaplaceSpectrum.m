%%
% Display of Laplacian spectrum on 2D domains.

name = 'france';
name = 'cat';
name = 'elephant';

names = {'cat'  'elephant'};

addpath('../toolbox/');
rep = MkResRep(name);


n = 100*2;
n0 = 400;
nc = 1000; % discretization (resampling) for the curves

% resample a curve using a fixed number of points
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );


Clist = {};
for it=1:2
    f = load_image(names{it}, n0);
    f = sum(f,3); f = double((f/max(f(:)))>.5);
    if f(1)==0
        f = 1-f;
    end
    [C,h] = contour(f,[.5,.5]);
    m = C(2,1);  c = C(:,2:m+1); c = (c-1)/(n0-1);
    Clist{end+1} = resample(c(1,:)'+1i*c(2,:)', nc);
end

% rotate the curves to match them
err = Inf; c = [];
for k=1:nc
    e = norm( Clist{1} - circshift(Clist{2},k) );
    if e<err
        err = e;
        c = circshift(Clist{2},k);
    end
    e = norm( Clist{1} - circshift(Clist{2}(end:-1:1),k) );
    if e<err
        err = e;
        c = Clist{2}(end:-1:1);
    end
end
Clist{2} = c;


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


q = 50; % ?frames

klist = [2];
klist = [2 3 4 5 10];
kmax = max(klist); % #eigenvectors
        
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
    
    
% display
r = 15; % #levellines
t = linspace(0,1,n);

Uold = [];
for it=1:q
    s = (it-1)/(q-1);
    % load shape
    c = Clist{1}*(1-s) + Clist{2}*s;
    %
    Img = 1 - inpolygon(X,Y,imag(c), real(c));
    % laplacian
    rho = Img+1e-5;
    Delta = Grad' * spdiags( [rho(:); rho(:)], 0, 2*n*n,2*n*n ) * Grad;
    % do SVD
    [U,S,V] = eigs(Delta,kmax+1, 'SM'); S = diag(S);
    [~,I] = sort(S); U = U(:,I); U = U(:,2:end);
    if it>1
        Z = diag(sign(diag(U'*Uold)));
        U = U*Z;
    end
    Uold = U;
    for i=klist
        u = reshape(U(:,i), [n,n]);
        u = u/max(abs(u(:)));
        
        % u(Img==1) = NaN;
        
        opt.cm = jet(r-1);
        B = display_shape_function( (u+1)/2,opt);
        I = find(Img==1);
        B([I,I+n*n,I+2*n*n])=1;
        
        
        % display
        clf; hold on;
        imagesc(t,t,B);
        contour(t,t,u,linspace(-1,1,r), 'k');
        caxis([-1 1]);
        plot(c([1:end 1]), 'k', 'LineWidth', 3);
        axis image; axis off; axis ij;
        drawnow;
        saveas(gcf, [rep 'eig-' znum2str(i,2) '-' znum2str(it,2) '.png' ]);
    end
end

 % AutoCrop(rep, 'eig');

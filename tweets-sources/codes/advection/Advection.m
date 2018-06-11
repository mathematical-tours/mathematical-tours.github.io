%%
% Advection PDE

n = 512; 

addpath('../toolbox/');
rep = MkResRep();

% gaussian blur
randn('state', 1234);
t = [0:n/2,-n/2+1:-1]';
G = @(s)exp(-t.^2/(2*s^2)); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));
normalize = @(V)V ./ repmat( max(1e-9,sqrt(sum(V.^2, 3))) , [1 1 2]);

% Laplacian
a = [2:n n]; b = [1 1:n-1]; % sym
a = [2:n 1]; b = [n 1:n-1]; % per
% grad
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f)/2;
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) )/2;
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

% image
name = 'hibiscus';
name = 'grid';
switch name
    case 'grid'
        [Y,X] = meshgrid(1:n,1:n);
        r = 32; k = 6;
        f0 = 1 - double( mod(X,r)<=k | mod(Y,r)<=k );
    otherwise
        f0 = rescale(sum(load_image(name, n),3));
        [fS,I] = sort(f0(:)); f0(I) = linspace(0,1,length(f0(:)));
end


% advection field
s = 80;
u = cat(3, GFilt(randn(n),s), GFilt(randn(n),s) );
%s = 80;
%u = Grad(  GFilt(randn(n),s) );
u = u/max(abs(u(:)));
% u = normalize(u);



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

% example of advection along the flow
rho = 3;
q = 70; %#display frames
k = 1;
f = f0;
for i=1:q
    clf; imagesc(f); colormap gray(256);
    axis image; axis off; axis xy;
    drawnow;
    % saveas(gcf, [rep 'evol-' znum2str(k,3) '.png'], 'png');
    imwrite(rescale(f), [rep name '-' znum2str(i,3) '.png']);
    f = W(f, rho*u );
end


% to display VF
sf = 1.2; % stretch factor
lw = 2;
q = 20;
k = round(n/q); % subsampling
plotvf = @(v,col)quiver(v(1:k:end,1:k:end,2),v(1:k:end,1:k:end,1), sf, 'Color', col, 'LineWidth', lw);

clf; hold on; 
        plotvf(u, 'b');
        axis equal; axis off;
        saveas(gcf, [rep 'flow.png'], 'png');


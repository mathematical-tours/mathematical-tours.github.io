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

% image
f0 = randn(n);


% advection field
s = 80;
s = 60;
u = cat(3, GFilt(randn(n),s), GFilt(randn(n),s) );
%s = 80;
%u = Grad(  GFilt(randn(n),s) );
u = u/max(abs(u(:)));
 u = normalize(u);



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
rho = 1/4;
q = 50; %#display frames
sk = 1;
k = .5;
f = f0; % equalize(f0);
g = f0;
F = f;
itsvg = 1;
for it=1:q*sk
    
    if 1 % mod(it,sk)==1
        F1 = equalize(F,f0);
        % F1 = F;
        
        clf; imagesc(F1); colormap gray(256);
        axis image; axis off; axis xy;
        caxis([-1 1]*3);
        drawnow;
        imwrite(rescale(clamp(F1,-3,3)), [rep  'anim-' znum2str(itsvg,3) '.png']);
        itsvg = itsvg+1;
    end
    
    
    %f = double(f>median(f(:)));
    
    
    % saveas(gcf, [rep 'evol-' znum2str(k,3) '.png'], 'png');
    f = W(f, rho*u );
    g = W(g,-rho*u );
    f = equalize(f,f0);
    g = equalize(g,f0);
    F = (1-1/it)*F + 1/it*(f+g)/2;
end


% to display VF
sf = 1.2; % stretch factor
lw = 2;
q = 30;
k = round(n/q); % subsampling
% plotvf = @(v,col)quiver(v(1:k:end,1:k:end,2),v(1:k:end,1:k:end,1), 'filled' 'Color', col, 'LineWidth', lw, 'AutoScaleFactor', 1.2);
clf; hold on; 
%plotvf(u, 'k');
quiver(u(1:k:end,1:k:end,2),u(1:k:end,1:k:end,1), 'k', 'filled', 'LineWidth', 1, 'AutoScaleFactor', .8);    
axis equal; axis off;
saveas(gcf, [rep 'flow.png'], 'png');


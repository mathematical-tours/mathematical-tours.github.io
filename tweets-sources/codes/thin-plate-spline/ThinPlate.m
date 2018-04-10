%%
% Test for thin plate spline interpolation.

addpath('../toolbox/');
rep = MkResRep();

%%
% Load images.

n = 256;
names = {'van-gogh' 'van-gogh-3'};
for k=1:2
    f{k} = rescale( load_image(names{k}, n) );
end

% pick points
C = distinguishable_colors(50);
clf; 
for k=1:2
    subplot(1,2,k);  hold on;
    imagesc(f{k}); axis image; axis off; axis ij;
end
A = []; B = []; k = 0;
while true
    k = k+1;
    %%%% 
    subplot(1,2,1); 
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 25, 'Color', C(k,:));
    if button==3
        break;
    end
    A(end+1) = a+1i*b;
    %%%%%
    subplot(1,2,2); 
    plot(a,b, 'x', 'MarkerSize', 20, 'Color', C(k,:));
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 25, 'Color', C(k,:));
    B(end+1) = a+1i*b;
end
% add the 4 corners
if 0
    A(end+1:end+4) = [1, n, 1i*n, n+1i*n];
    B(end+1:end+4) = [1, n, 1i*n, n+1i*n];
end
    
A = A(:); B = B(:);
q = length(A); % #points

% pairwise distance matrix
A1 = [real(A) imag(A) ones(q,1)]'; % use homogeneous coords
D = distmat(A1);
% kernel function
phi = @(u)u.^2 .* log(u+1e-10);
K = phi(D);
% linear system 
L = [K, A1'; A1, zeros(3,3)];
b = [real(B) imag(B)]; b(end+1:end+3,:) = 0;
S = L\b;
% weights
W = S(1:end-3,:);
% affine coefs
Wx = S(end-2,:);
Wy = S(end-1,:);
Wc = S(end-0,:);

% warp the position
Kern = @(X,Y)distmat([X(:) Y(:)]', [real(A) imag(A)]');
Extrap = @(X,Y) phi(Kern(X,Y)) * W + X(:)*Wx + Y(:)*Wy + ones(length(X),1)*Wc;

% generate a coarse horizontal grid
m0 = 256; m1 = 16;
%
[X1,Y1] = meshgrid(linspace(1,n,m0),linspace(1,n,m1));
XY1 = Extrap(X1(:),Y1(:));
XY1 = reshape(XY1, [m1 m0 2]);
%
[X2,Y2] = meshgrid(linspace(1,n,m1),linspace(1,n,m0));
XY2 = Extrap(X2(:),Y2(:));
XY2 = reshape(XY2, [m0 m1 2]);
% input shape
clf; hold on; 
imagesc(.3 + f{1}); axis image; axis off;
for k=1:m1
    plot(X1(k,:), Y1(k,:), 'k', 'LineWidth', 2);
    plot(X2(:,k), Y2(:,k), 'k', 'LineWidth', 2);
end
plot(real(A), imag(A), 'b.', 'MarkerSize', 30);
axis ij;
saveas(gcf, [rep 'input.png']);
% Target shape
clf; hold on; 
imagesc(.3 + f{2}); axis image; axis off;
for k=1:m1
    plot(XY1(k,:,1), XY1(k,:,2), 'k', 'LineWidth', 2);
    plot(XY2(:,k,1), XY2(:,k,2), 'k', 'LineWidth', 2);
end
plot(real(B), imag(B), 'r.', 'MarkerSize', 30);
axis ij;
saveas(gcf, [rep 'output.png']);



% interpolate at grid position
[X,Y] = meshgrid(1:n,1:n);
XY1 = Extrap(X(:),Y(:));
F = [];
for k=1:3
    u = interp2(X,Y,f{2}(:,:,k),XY1(:,1),XY1(:,2));
    F(:,:,k) = reshape(u, [n n]);
end

clf;
imagesc(F); axis image; axis off;


r = 40;
for k=1:r
    t = (k-1)/(r-1);
    imwrite(rescale( (1-t)*f{1}+t*F ), [rep 'interp-' znum2str(k,2) '.png']);
end

% linear blending
r = 30;
while true
for k=[1:r, r-1:-1:2]
    t = (k-1)/(r-1);
    clf; imagesc((1-t)*f{1}+t*F); axis image; axis off;
    drawnow;
end
end


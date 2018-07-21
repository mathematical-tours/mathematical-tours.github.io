%%
% Edge detection classical technics.

addpath('../toolbox/');
rep = MkResRep();

% image size
n = 256*2;

name = 'hibiscus';
f = load_image(name, n);
f = rescale(sum(f,3));
[fS,I] = sort(f(:)); f(I) = linspace(0,1,length(f(:)));


% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

% operators
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f)/2;
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) )/2;
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

%
dx = @(f)(f(a,:)-f(b,:)) / 2;
dy = @(f)dx(f')';
d2x = @(f)f(a,:) + f(b,:) - 2*f;
d2y = @(f)d2x(f')';
dxy = @(f)dy(dx(f));
Hessian = @(f)cat(3, d2x(f), dxy(f), d2y(f));

% imwrite(rescale(f), [rep 'original.png']);


% display quantized colormap
t = linspace(0,1,n);
r = 12; % #levellines


q = 70; 
slist = linspace(1,15, q);
for i=1:q
    s = slist(i);
    fs = GFilt(f,s);
    %
    if 0
    clf; hold on;
    imagesc(t,t,fs);
    contour(t,t,fs,linspace(0,1,r), 'k');
    colormap(gray(r-1));
    caxis([0 1]);
    axis image; axis off;
    drawnow;
    end
    
    imwrite(fs, [rep 'smoothed-' znum2str(i,2) '.png']);
    
    % laplacian levelset
    D = Delta(fs);
    clf; hold on;
    imagesc(t,t,rescale(clamp(rescale(D),.2,.8)));
    contour(t,t,D,[0 0], 'b', 'LineWidth', 2);
    colormap(gray(256));
    caxis([0 1]);
    axis image; axis off;
    % drawnow;
    saveas(gcf,[rep 'laplacian-' znum2str(i,2) '.png']);

    
    % Canny
    g = Grad(fs);
    h = Hessian(fs);
    a = h(:,:,1:2).*repmat(g(:,:,1), [1 1 2]) + ...
        h(:,:,2:3).*repmat(g(:,:,2), [1 1 2]);
    % < H*(nabla fs),nabla fs>
    D = sum(a.*g, 3);    
    %  levelset
    clf; hold on;
    imagesc(t,t,rescale(clamp(rescale(D),.2,.8)));
    contour(t,t,D,[0 0], 'r', 'LineWidth', 2);
    colormap(gray(256));
    caxis([0 1]);
    axis image; axis off;
	saveas(gcf,[rep 'canny-' znum2str(i,2) '.png']);
    drawnow;
end

AutoCrop(rep, ['laplacian-']);

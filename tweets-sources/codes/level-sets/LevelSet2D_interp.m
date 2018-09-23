%%
% Display of levelsets of the distance function


addpath('../toolbox/');
rep = MkResRep('interp2d');

names = {'elephant' 'cat'};
names = {'bunny' 'elephant'};


n = 200; % size of the image for loading


% gaussian blur
t = [0:n/2,-n/2+1:-1]';
G = @(s)exp(-t.^2/(2*s^2)); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));



% 
for s=1:2
f = load_image(names{s}, n);
f = rescale(sum(f,3))>.5;
if f(1)==1
    f = 1-f;
end


% extract a level set as a planar curve
f1 = GFilt(f,1.5);
[C,h] = contour(f1,[.5,.5]);
m = C(2,1);  c{s} = C(:,2:m+1);

% distance transform
[X,Y] = meshgrid(1:n,1:n);
D = distmat(c{s}, [X(:) Y(:)]');
[D,S] = min(D,[],1);
D = reshape(D,[n n]);
S = reshape(S,[n n]);
D(f1<.5) = -D(f1<.5);
Dlist{s} = D;
end

q = 50;

for i=1:q
    t = (i-1)/(q-1);
    D = Dlist{1}*(1-t) + Dlist{2}*t;
    U = clamp(D/max(D(:)),-1,1);
    r = 17;
    clf; hold on;
    imagesc(1:n,1:n,U);
    % plot(c(1,:), c(2,:), 'r', 'LineWidth', 2);
    contour(1:n,1:n,U,linspace(-1,1,r), 'k');
    contour(1:n,1:n,U,[0 0], 'k', 'LineWidth', 3);
    axis equal; axis([1 n 1 n]); axis off; axis ij;
    colormap(jet(r-1));
    caxis([-1 1]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, 'anim-');
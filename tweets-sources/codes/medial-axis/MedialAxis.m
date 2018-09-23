%%
% Compute distance and medial axis.

addpath('../toolbox/');
rep = MkResRep();

name = 'elephant';
name = 'cat';


n = 400; % size of the image for loading

% 
f = load_image(name, n);
f = rescale(sum(f,3))>.5;
if f(1)==1
    f = 1-f;
end


% extract a level set as a planar curve
[C,h] = contour(f,[.5,.5]);
m = C(2,1);  c = C(:,2:m+1);

% distance transform
[X,Y] = meshgrid(1:n,1:n);
D = distmat(c, [X(:) Y(:)]');
[D,S] = min(D,[],1);
D = reshape(D,[n n]);
S = reshape(S,[n n]);
D(f==0) = 0;

% convert to color image
CM = parula(256);
I = 1+floor(rescale(D)*255);
U = reshape(CM(I,:), [n n 3]);
for i=1:3
    u = U(:,:,i); u(f==0) = 1; U(:,:,i) = u;
end

clf; hold on;
imagesc(1:n,1:n,U);
plot(c(1,:), c(2,:), 'r', 'LineWidth', 2);
contour(1:n,1:n,D,linspace(0,max(D(:)),10), 'k');
axis equal; axis([1 n 1 n]); axis off; axis ij;
saveas(gcf, [rep name '-dist.png']);

% compute gradient
q = size(c,2);
i = [1 1:n-1]; j = [2:n n];
mymod = @(s) s.*( s<q/2 & s>-q/2 ) + (s-q).*(s>=q/2)  + (s+q).*(s<=-q/2);
Dx = @(f)mymod(f(j,:) - f(i,:));
Dy = @(f)Dx(f')';
E = sqrt(Dx(S).^2+Dy(S).^2);
E(f==0) = 0;


t = 20;
clf; hold on;
imagesc(1:n,1:n,1-(E>t));
plot(c(1,:), c(2,:), 'r', 'LineWidth', 2);
axis equal; axis([1 n 1 n]); axis off; axis ij;
saveas(gcf, [rep name '-skel.png']);








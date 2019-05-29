%%
% Test for the Durand?Kerner method.

if not(exist('test'))
    test=1;
end

addpath('../toolbox/');
rep = MkResRep(num2str(test));

% n = 5;
% z0 = randn(n,1)+1i*randn(n,1);

% click and play
clf; hold on;
z0 = [];
while true
    axis equal; axis([-1 1 -1 1]); 
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    if button==3
        break;
    end
    z0(end+1) = a+1i*b;
end
z0 = z0(:);
n = length(z0);
% click and play
clf; hold on;
plot(z0, 'k.', 'MarkerSize', 25);
z1 = [];
for i=1:n
    axis equal; axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    z1(end+1) = a+1i*b;
end
z1 = z1(:);


niter = 600;
niter = 150*4;

tau = .02;

% z = randn(n,1)+1i*randn(n,1);
Z = [];
z = z1;
for i=1:niter
    Z(:,end+1) = z;
    q = z;
    for k=1:n
        J = setdiff(1:n,k);
        q(k) = z(k) - polyval(z(k),z0) ./ polyval(z(k),z(J));
    end
    z = (1-tau)*z + tau*q;
    
    D = abs( repmat(z,[1 length(z0)]) - repmat(permute(z0,[2 1]),[length(z) 1]) );
	d = mean(min(D));
    if d<=1e-2
        Z(:,end+1) = z;
        break;
    end
end

% evaluate polynomial on a grid
m = 256;
t = linspace(-1,1,m);
[Y,X] = meshgrid(t,t); XY = X+1i*Y;
F = ones(m);
for i=1:length(z0)
    F = F .* (XY-z0(i));
end
U = log(.001+abs(F));

% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(t,t,U');
contour(t,t,U',linspace(min(U(:)),max(U(:)),r), 'k');
colormap(parula(r-1));
caxis([min(U(:)) max(U(:))]);
axis image; axis off;
plot(z0, 'r.', 'MarkerSize', 25);
%
q = 50;
ndisp = round(linspace(1,size(Z,2)-1,q)); kdisp = 1;
for i=1:size(Z,2)-1
    s = (i-1)/(size(Z,2)-2);
    plot( permute(Z(:,i:i+1), [2 1]), '-', 'Color', [s 0 1-s], 'LineWidth', 3 );
    if ndisp(kdisp)==i
        axis tight; axis equal; axis([-1 1 -1 1]); 
        axis off;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(kdisp,2) '.png']);
        kdisp = kdisp + 1;
    end
end

% AutoCrop(rep, 'anim')

% test = test+1;
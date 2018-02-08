%%
% Contrast iterative mean vs iterative median.

rep = '../results/mean-median/';
[~,~] = mkdir(rep);

rescale = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));

% number of quantized values
m = 6;

% # nn
w = 3;
w1 = 2*w+1;
% size of image
n = 256;

[Y,X] = meshgrid(1:n,1:n);
[dY,dX] = meshgrid(-w:w,-w:w);
dX = reshape(dX, [1 1 w1 w1]);
dY = reshape(dY, [1 1 w1 w1]);
X = repmat(X, [1 1 w1 w1]) + repmat(dX, [n n 1 1]);
Y = repmat(Y, [1 1 w1 w1]) + repmat(dY, [n n 1 1]);
% bc
X = mod(X-1,n)+1;
Y = mod(Y-1,n)+1;
% patch extract
Pi = @(f)reshape( f(X + (Y-1)*n), [n n w1*w1] );

f0 = double(randn(n)>0);

f0 = floor(rand(n)*m)/(m-1);

name = 'mean';
name = 'median';

switch name
    case 'mean'
        oper = @(u)mean(u,3);
    case 'median'
        oper = @(u)median(u,3);
end

niter = 80;
f = f0;
for i=1:niter
    imagesc(f);     
    axis image; axis off; 
    colormap parula(256);
    % caxis([0,1])
    drawnow;
    imwrite(rescale(f), [rep name '-' num2str(i) '.png']);
    f  = oper(Pi(f));
end

return;

niter = 500;
nbdisp = 50;
ndisp = round( linspace(1,sqrt(niter),nbdisp).^2 );
ndisp = 1:niter;

q = max(abs(f0));
f = f0;
a = 1;
for it=1:max(ndisp);
    f = oper(f(I));
    if it<=length(ndisp) && it==ndisp(a)
        c = (a-1)/(nbdisp-1);
        clf;
        plot(f, 'LineWidth', 2, 'Color', [c 0 1-c]);
        axis([1 n -q q]); axis off;
        drawnow;
        a = a+1;
    end
end
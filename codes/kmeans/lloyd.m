%%
% test for k-means on a dense lattice.

name = 'gaussian';
name = 'linear';
name = 'uniform';
name = 'half';

rep = ['./results/lloyd/' name '/'];
[~,~] = mkdir(rep);

n = 400; % size of the image
N = n*n;

t = linspace(0,1,n);

[v,u] = meshgrid(t,t);
X = [u(:)';v(:)'];

% measure on the image
switch name
    case 'uniform'
        mu = ones(N,1); mu(1)=1+1e-3;
    case 'gaussian'
        m = [.8 .2];
        s = .15;
        mu = .01 + exp( - ( (X(1,:)-m(1)).^2 + (X(2,:)-m(2)).^2 )/(2*s^2) );
    case 'half'
        mu = 1+(X(1,:)<.5)*10; mu(1)=max(mu)+5;
    case 'linear'
        mu = 1-X(1,:);
end
mu = mu(:);


if not(exist('Y'))
    % click and play
    Y0 = [];
    clf; hold on;
    while true
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 20);
        if button==3
            break;
        end
        Y0(:,end+1) = [a;b];
    end    
end

% random points
k = 100;
Y0 = randn(2,k)*.08 + .8;

k = size(Y0,2);


niter = 500;
disp_list = [1:10 15:5:100 150:50:niter];
kdisp = 1;

Y = Y0;
for it=1:niter
    % voronoi cells
    D = distmat(X,Y);
    [D1,I] = min(D,[], 2);
    % display
    if it==disp_list(kdisp)
        kdisp = kdisp+1;
        J = reshape(I, [n n]);
        D1 = reshape(D1, [n n]);
        clf; hold on;
        imagesc(t,t,-reshape(mu, [n n])); colormap gray(256);        
        for l=1:k
            contour(t,t,J==l,[.5 .5], 'k', 'LineWidth', 2);
        end
        caxis([min(-mu(:)) max(-mu(:))])
        plot(Y(2,:), Y(1,:), 'r.', 'MarkerSize', 20);
        plot([0 1 1 0 0], [0 0 1 1 0], 'k', 'LineWidth', 2);
        axis off;  drawnow;
        saveas(gcf, [rep 'iter-' num2str(it)], 'png');
    end
    % New centroids
    for l=1:k
        for r=1:2
            Y(r,l) = sum(X(r,I==l) .* mu(I==l)') / sum(mu(I==l));
        end
    end
end
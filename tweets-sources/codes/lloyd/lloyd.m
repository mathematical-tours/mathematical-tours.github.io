%%
% test for k-means on a dense lattice.

addpath('../toolbox/');
rep = MkResRep();

name = 'gaussian';
name = 'linear';
name = 'half';
name = 'uniform';


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


if not(exist('Y0'))
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
% k = 100;
% Y0 = randn(2,k)*.08 + .8;

k = size(Y0,2);


repsvg = [rep name '-' num2str(k) '/'];
[~,~] = mkdir(repsvg);

niter = 500;
disp_list = [1:10 15:5:100 150:50:niter];
kdisp = 1;

% coloring
C = distinguishable_colors(k+1);
C(4,:) = [];

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
        if strcmp(name, 'uniform')
            % render in colors
            R = zeros(n,n,3);
            for m=1:k
                for ic=1:3
                    R(:,:,ic) = R(:,:,ic) + (J==m)*C(m,ic);
                end
            end
            s = .5;
            imagesc(t,t,s + (1-s)*R); colormap jet(256);
        else
            % display density
            imagesc(t,t,-reshape(mu, [n n])); colormap gray(256); 
            caxis([min(-mu(:)) max(-mu(:))]);
        end
        for l=1:k
            contour(t,t,J==l,[.5 .5], 'k', 'LineWidth', 2);
        end
        for m=1:k
            plot(Y(2,m), Y(1,m), '.', 'MarkerSize', 25, 'Color', .8*C(m,:));
        end
        axis equal; 
        plot([0 1 1 0 0], [0 0 1 1 0], 'k', 'LineWidth', 2);
        axis off;  drawnow;
        saveas(gcf, [repsvg 'iter-' znum2str(it,3)], 'png');
    end
    % New centroids
    for l=1:k
        for r=1:2
            Y(r,l) = sum(X(r,I==l) .* mu(I==l)') / sum(mu(I==l));
        end
    end
end

% AutoCrop(repsvg, 'iter-') 
% convert iter-*.png lloyd.gif 
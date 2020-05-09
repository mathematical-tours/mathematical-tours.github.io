%%
% Display of K-NN classification method.

addpath('../toolbox/');
rep = MkResRep();

% generate points which will be center of
if not(exist('X'))
    % click and play
    X = {}; col = 'b';
    clf; hold on;
    for i=1:2
        X{i} = [];
        axis equal;
        while true
            axis([0 1 0 1]);
            [a,b,button] = ginput(1);
            plot(a,b, [col '.'], 'MarkerSize', 20);
            if button==3
                break;
            end
            X{i}(:,end+1) = [a;b];
        end
        col = 'r';
    end
end

% re-sample with noise
s = .06;
n = 200; % #sample per distributions
for i=1:2
    p = size(X{i},2);
    Y{i} = [];
    for j=1:n
        r = ceil(rand*p);
        Y{i}(:,end+1) = X{i}(:,r) + randn(2,1)*s;
    end
end


% Knn on a grid
q = 256*2;
u = linspace(0,1,q);
[V,U] = meshgrid(u,u);
Z = [U(:) V(:)]';
D = distmat(Z,[Y{1},Y{2}]);
[~,I] = sort(D,2);
Ic = (I<=n);
%
NNk = @(k)reshape(sum(Ic(:,1:k),2)/k, [q q]);

col = {'b' 'r'};
ms = 15;
klist = 1:50;
for i=1:length(klist)
    k = klist(i);

    R = NNk(k);
    %
    T = (1/2:1:k-1/2)/k;
    if k==1
        T = [1/2 1/2];
    end
    %
    m = linspace(0,1,k+1)';
    CM = m*[.5 .5 1] + (1-m)*[1 .5 .5];
    %
    clf; hold on;
    imagesc(u,u,R');
    contour(u,u,R', T ,'k'); %
    colormap(CM);
    caxis([0 1]);
    axis image; axis off;
    %
    for i=1:2
        plot(Y{i}(1,:), Y{i}(2,:), '.', 'color', col{i}, 'MarkerSize', ms);
    end
    axis equal; axis([0 1 0 1]); axis off;
    drawnow;
	saveas(gcf, [rep 'anim-' znum2str(k,2) '.png'], 'png');
end

% AutoCrop(rep, 'anim-'); 


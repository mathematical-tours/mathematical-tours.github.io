%%
% Doing displacement interpolation for color images.


addpath('../toolbox/');
addpath('./toolbox-lsap/');

names = {'cezanne' 'multi-3'};
name = [names{1} '-' names{2}];

rep = ['../results/color-displacement/' name '/'];
[~,~] = mkdir(rep);

n = 100; % low res
N = 240; % high res
f = {}; F = {}; X = {}; 
for i=1:2
    f{i} = rescale( load_image(names{i}, n) );
    F{i} = rescale( load_image(names{i}, N) );
    X{i} = reshape(f{i}, [n*n 3]);
end

% compute OT
C = distmat(X{1}',X{2}').^2;
C1 = int32( round(C*1e6) );
tic;
[J,varrho,u,v] = hungarianLSAP(C1);
toc
% Reverse permutation
I = ones(n*n,1); I(J) = 1:n*n;

m = 50; 
for it=1:m
    t = (it-1)/(m-1);
    col = [1-t;0;t];    
    Xt{1} = (1-t)*X{1}      + t*X{2}(J,:);
    Xt{2} = (1-t)*X{1}(I,:) + t*X{2}; % same cloud but other indexing 
    % display cloud    
    s = ones(n*n,1)*20; % size
    clf;
    scatter3(...
        Xt{1}(:,1), ...
        Xt{1}(:,2), ...
        Xt{1}(:,3),s, Xt{1}, 'filled' );
    axis ij; view(3); drawnow;
    axis equal; axis([0 1 0 1 0 1]);
    box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    drawnow;
    saveas(gcf, [rep 'interp-' znum2str(it,2) '.png']);
    % display image
    for k=1:2
        ft = reshape(Xt{k}, [n n 3]);
        Ft = image_resize(ft,[N N]) + (F{k}-image_resize(f{k}, [N N])); % add back details
        % upscale    
        imwrite(clamp(Ft), [rep names{k} '-' znum2str(it,2) '.png'], 'png');
    end
end

% ask for cropping
if 0
    AutoCrop(rep, 'interp-')
    % convert interp-*.png interp.gif
    % convert cezanne-*.png cezanne.gif
    % convert multi-3-*.png multi-3.gif
end
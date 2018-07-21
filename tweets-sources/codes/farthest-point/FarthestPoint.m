%%
% Test for Farthest point sampling.


addpath('../toolbox/');
rep = MkResRep();

n = 256;
name = 'hibiscus';
f = rescale(sum(load_image(name, n),3));
% equalize
[fS,I] = sort(f(:)); f(I) = linspace(0,1,length(f(:)));

% points
P = [1;1];

[Y,X] = meshgrid(1:n,1:n);
Q = [X(:)';Y(:)'];

fmin = .05;
U = rescale(-f,0.01,1).^2;

pmax = 400;

q = 50; % number of displayed
ndisp = round(linspace(1,pmax,q));
idisp = 1;

Col = distinguishable_colors(pmax);


for i=1:pmax
    %
    D = distmat(P,Q);
    if 0
        w = U(P(1,:) + (n-1)*P(2,:))+fmin;
        D = D .* repmat(w(:), [1 n*n]);
    else
        D = D .* repmat(U(:)', [size(D,1) 1]);
    end
    [D,Qi] = min(D,[],1);
    Qi = reshape(Qi,[n n]);
    [~,j] = max(D(:));
    %
    if ndisp(idisp)==i
        % generate color back
        C = zeros(n,n,3);
        for s=1:3
            for m=1:i
                C(:,:,s) = C(:,:,s) + Col(m,s) .* (Qi==m);
            end
        end
        % Modulate
        C = (.4+.6*C) .* (repmat(f, [1 1 3]));
        % display image with voronoi background
        clf; hold on;
        imagesc(1:n,1:n,C);
        s = ones(i,1)*30; % size
        scatter( P(2,:), P(1,:), s, Col(1:i,:), 'filled' );
        axis image; axis off;
        drawnow;
        saveas(gcf, [rep 'sample-' znum2str(idisp,2) '.png']);
        idisp = idisp+1;
    end
    
    P(:,end+1) = [X(j);Y(j)];
end

% AutoCrop(rep, 'sample-');
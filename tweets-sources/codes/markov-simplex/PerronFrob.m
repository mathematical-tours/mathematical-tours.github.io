%%
% Test for contraction for stochastic matrices (aka Perron-Frobenius
% theomr).

addpath('../toolbox/');
rep = MkResRep();

% dimension
n = 3; 

u = ones(n,1);

name = 'test2';
name = 'test2';
name = 'isotropic';
name = 'test3';
name = 'rand';

% input matrix
switch name 
    case 'rand'
        K = rand(n)+20*eye(n);
    case 'isotropic'
        K = ones(n)+20*eye(n);
    case 'test1'
        K = [20 1 4; 5 19 1; 3 1 21];
    case 'test2'
        K = [16 5 4; 5 20 0; 4 0 21];
    case 'test3'
        K = [0 0 1; 0 0 1; 1 0 0];
end

tau = .9;
K = tau*eye(3)+(1-tau)*K;

% make it stochastic, K'*u=1
K = K*diag(1./(K'*u));


% display tiangle
a = exp(1i*pi/2 + 2i*pi*(0:n-1)'/n);

% invariant proba
[V,S] = eig(K); 
[~,I] = sort(diag(S)); 
v = V(:,I(end)); v=v/sum(v);

lw = 2; ms = 30;

q = 40;
for i=1:q
    X = eye(n); % each column is a point of the simplex
    clf; hold on;
    for j=1:i
        t = (j-1)/(q-1);
        z = X'*a; % barycentric coordinates
        h = plot(z([1:end 1]), '-', 'MarkerSize', ms, 'LineWidth', lw, 'color', [t 0 1-t]);
        if j<i
            h.Color(4) = .15; % set(h,'alpha', .5);
        end
        X = K*X;
    end
    plot(sum(v.*a), 'k.', 'MarkerSize', ms);
    axis tight;  axis equal; axis off;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,2) '.png']);
end

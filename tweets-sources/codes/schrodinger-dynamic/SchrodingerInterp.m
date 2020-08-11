%%
% test to display Brownian bridge interpolation.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

randn('state', 66);
% click selection
if not(exist('X'))
X = []; Y = [];
clf; hold on;
while true
    axis equal; axis([0 1 0 1]); box on; set(gca, 'XTick', [], 'YTick', []);
    [a,b,button] = ginput(1);
    if button==1
        X(:,end+1) = [a;b];
        plot(a,b, 'b.', 'MarkerSize', 15);
    elseif button==3
        Y(:,end+1) = [a;b];
        plot(a,b, 'r.', 'MarkerSize', 15);        
    else
        break;
    end
end
end
N = size(X,2);
M = size(Y,2);

P = 256*2; % number of step for brownian bridge simulation
K = 50; % numbef of bridges

% 
c = distmat(X,Y).^2;
mu = ones(N,1)/N; 
nu = ones(M,1)/M; 
% budget of path
Np = 100; 

epsilon = .02;
epsilon = .001;
epsilon = .07;
epsilon = .11;
options.niter = 20000;
[u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);

%%
% Display results using Brownian Bridges

sigma = .08*sqrt(epsilon);
sigma = .2*sqrt(epsilon);
if epsilon<=1e-3
    sigma = 0;
end

G = round(gamma*Np);
[I,J,NP] = find(G);

Bc = [];
for s=1:length(I)
    a = X(1,I(s)) + 1i * X(2,I(s));
    b = Y(1,J(s)) + 1i * Y(2,J(s));
    K = NP(s);
    % generate a bunch of bridges
    Bc(:,end+1:end+K) = brownian_bridge(P,K,a,b,sigma);
end


q = 50;
for it=1:q
    s = (it-1)/(q-1);
    I = ceil(1+s*(P-1));
    clf; hold on;
    l = 10;
    for j=1:l:I-1
        t = (j-1)/(P-1);
        plot(real(Bc(j:min(j+l,I),:)),imag(Bc(j:min(j+l,I),:)),'LineWidth', 1, 'color', .5*[1-t 0 t]+.5);
    end
    plot(real(Bc(I,:)),imag(Bc(I,:)), '.', 'color', [1-s 0 s], 'MarkerSize', 15);
    plot(X(1,:), X(2,:), 'r.', 'MarkerSize', 25);
    plot(Y(1,:), Y(2,:), 'b.', 'MarkerSize', 25);
    axis off; axis equal;
    axis([0 1 0 1]);
    drawnow;
    saveas(gcf, [rep 'schrodinger-' znum2str(it,2) '.png'], 'png');    
end
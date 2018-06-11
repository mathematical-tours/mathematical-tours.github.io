%%
% test to display Brownian bridge interpolation.


addpath('../toolbox/');
rep = MkResRep();

% number of points
N = 6;


randn('state', 66);
X = randn(2,N);
Y = randn(2,N)+.5;
% click to get points

% click selection
X = []; Y = [];
clf; hold on;
while true
    axis([0 1 0 1]);
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
N = size(X,2);
M = size(Y,2);



P = 256; % number of step for brownian bridge simulation
K = 50; % numbef of bridges


a = 0+0*1i;
b = 1+1*1i;
sigma = .05;
Bc = brownian_bridge(P,K,a,b,sigma);

clf;
plot_colored(Bc);

%%
% Run Sinkhorn


c = distmat(X,Y).^2;
mu = ones(N,1)/N; 
nu = ones(M,1)/M; 


Np = 100; % budget of path

q = 50;
elist = linspace(.005, .2, q);
elist = .003 + .4*linspace(0,1,q).^2;

for it=1:q
    epsilon = elist(it);
    
    options.niter = 5000;
    [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);
    
    %%
    % Display results using Brownian Bridges
    
    sigma = .08*sqrt(epsilon);
    if epsilon<.01
        sigma=0;
    end
    
    G = round(gamma*Np);
    [I,J,NP] = find(G);
    clf; hold on;
    for s=1:length(I)
        a = X(1,I(s)) + 1i * X(2,I(s));
        b = Y(1,J(s)) + 1i * Y(2,J(s));
        K = NP(s);
        % generate a bunch of bridges
        Bc = brownian_bridge(P,K,a,b,sigma);
        % plot
        plot_colored(Bc);
    end
    plot(X(1,:), X(2,:), 'r.', 'MarkerSize', 25);
    plot(Y(1,:), Y(2,:), 'b.', 'MarkerSize', 25);
    axis off; axis equal;
    axis([0 1 0 1]);
    drawnow;
    saveas(gcf, [rep 'schrodinger-' znum2str(it,2) '.png'], 'png');    
end
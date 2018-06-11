
addpath('../toolbox/');

%% Bridge

rep = MkResRep('bridge');


P = 256; % number of step for brownian bridge simulation
K = 50; % numbef of bridges

a = (1+1i)*.1;
b = (1+1i)*.9;

q = 20;
slist = linspace(0,.2,q);

for it=1:q
    sigma = slist(it);
    Bc = brownian_bridge(P,K,a,b,sigma);
    clf;
    plot_colored(Bc);
    plot(a, 'r.', 'MarkerSize', 25);
    plot(b, 'b.', 'MarkerSize', 25);
    axis off; axis equal;
    axis([0 1 0 1]);
    drawnow;
end


%% Brownian

rep = MkResRep('bridge');

P = 256; % number of step for brownian bridge simulation
K = 50; % numbef of bridges


q = 20;
slist = linspace(0,.2,q);

for it=1:q
    sigma = slist(it);
    
    t = linspace(0,1,P)';
    % generate a bunch of bridges
    B = (randn(P-1,K)+1i*randn(P-1,K))/sqrt(P);
    B = [zeros(1,K); cumsum(B)];
    B = B*sigma;

    clf;
    plot_colored(B);
    plot(0,0, 'r.', 'MarkerSize', 25);
    axis off; axis equal;
    axis([-1 1 -1 1]);
    drawnow;
end
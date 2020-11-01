%%
% display of invariant measures.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);



% logistic map.
m = 2;  % single stable point
m = 4; % fully chaotic
m = 3.8;  % intermediate
m = 2.5;

% lagrangian discretization
ndiscr = 10000;
ndiscr = 5000000;
Z = linspace(0,1,ndiscr); % deterministic init
Z = rand(ndiscr,1); % random init
p = 200; vmax = 8;
t = linspace(0,1,p);

q = 75; % anim plots
niter = 60;
mlist = linspace(3.5,4,q);
for it=1:q
    m = mlist(it);
    f = @(x)m*x.*(1-x); 
    s = (it-1)/(q-1);    
    Z = linspace(0,1,ndiscr); % deterministic init
    for iter=1:niter
        Z = f(Z);
    end
	h = hist(Z,t);
	h = p*h/sum(h);
    clf; area(t, h, 'FaceColor', [s 0 1-s], 'EdgeColor', [s 0 1-s], 'LineWidth', 2);
    axis([0 1 0 vmax]);
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);    
end
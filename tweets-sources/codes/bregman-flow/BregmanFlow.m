%%
% Display of Bregman flow, i.e. flow on a hessian manifold.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


psi = @(s,a) (s.^a-a.*(s-1)-1) ./ (a.*(a-1));

x = linspace(0,4,200)';

alist = [-1 0 1 2] + 1e-5;
alist = (-1:.5:3)  + 1e-5;
alist = (-1:1:2)  + 1e-5;

alist = (-2:.5:3)  + 1e-5;

lgd = {};
clf; hold on;
for j=1:length(alist)
    a = alist(j);
    s = (j-1)/(length(alist)-1);
    plot( x, psi(x,a), 'color', [s 0 1-s], 'LineWidth', 2 );    
    lgd{end+1} = ['\alpha=' num2str(round(a*4)/4)];
end
axis([0 max(x) 0 5]);
set(gca, 'PlotBoxAspectRatio', [1 2/3 1])
box on;  set(gca,'XTick', [],'YTick', []);
saveas(gcf, [rep 'Entropies.eps'], 'epsc');
legend(lgd);

% in 2^d, |A*x-y|^2
A = [1.2,1];
y = 1; 
nablaf = @(x)A'*(A*x-y);
a = 0;

niter = 2000;
tau = .05; 

beta_list = 2-alist;

q = 200;
for it=1:q
    
    t = (it-1)/q;
        
    clf; hold on;
    for j=1:length(beta_list)
        beta = beta_list(j);
        s = (j-1)/(length(beta_list)-1);
       
        x = max(1e-4, .52 + .5 * [cos(2*pi*t); sin(2*pi*t)]);
        for jt=1:niter
            d = tau * (x(:,end).^beta) .* nablaf(x(:,end));
            if norm(d)>.05
                d = .05 * d/norm(d);
            end
            x(:,end+1) = x(:,end) - d;
        end
        plot([0 1], [1 0], ':', 'color', [0 .5 0], 'LineWidth', 2);
        plot(x(1,:),x(2,:), 'LineWidth', 2, 'color', [s 0 1-s])
        plot(x(1,1),x(2,1), 'k.', 'MarkerSize', 20);
        plot(x(1,end),x(2,end), '.', 'MarkerSize', 20, 'color', [s 0 1-s]);
        axis equal; axis([0 1.05 0 1.05]); box on;
        set(gca,'XTick', [],'YTick', []);
    end
    drawnow;
    % mysaveas(it);
end






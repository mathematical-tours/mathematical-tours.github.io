%%
% Mean shift method.

if not(exist('test'))
    test = 0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));

% Sample points.
s = .05*3; % spread
p = 30*20; % # point each time
x = [];
clf; hold on;
while true
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    [a,b,button] = ginput(1);
    if button==3
        break;
    end
    x = [x; a+1i*b + s*randn(p,1) + 1i*s*randn(p,1)];
    plot(x, 'k.', 'MarkerSize', 15);
end

% kernel
sigma = .08;
c2v = @(u)[real(u) imag(u)];
K = @(x,y)exp(-distmat(c2v(x)',c2v(y)').^2/(2*sigma^2));

niter = 200;
q = 50; 
ndisp = round(linspace(1,niter,q));
kdisp = 1;

y = x(1:1:end);
Y = [y];
tau = .04;
for it=1:niter
    t = (it-1)/(niter-1);
    % 
    y = (1-tau)*y + tau * (K(y,x)*x) ./ sum(K(y,x),2);
    Y(:,end+1) = y;
    if it==ndisp(kdisp)
        clf; hold on;
        plot(x, 'k.', 'MarkerSize', 15);
        plot(transpose(Y), '-', 'Color', [t 0 1-t]);
        plot(y, '.', 'Color', [t 0 1-t], 'MarkerSize', 30);
        axis equal; axis([0 1 0 1]); axis off;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(kdisp,2) '.png']);
        kdisp = kdisp+1;
    end
end

% AutoCrop(rep, 'anim');
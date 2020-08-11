%%
% Use a reweighting l2 scheme for the computation of medians.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 50;
p = round(.16*n); % outliers

randn('state', 1235);
crandn = @(n)randn(n,1) + 1i*randn(n,1);

x = [crandn(n-p)*3; ...
        crandn(p)*.3 + 30 *( .4+1i*.6 ); ...
        crandn(p)*.3 + 30 *( .6+1i*.6 ); ...
        crandn(p)*.3 + 30 *( .4+1i*.8 )];
y = x;


niter = 70; 




tau = .15; % slow down
p = 1;

tau = .08; % slow down
p = 2;

tau = .025; % slow down 
p = 6;

epsilon = 1e-9;
epsilon = 1e-5;
clf; hold on;
for it=1:niter
    s = (it-1)/(niter-1);
    ynew = [];
    for i=1:length(y)
        w = 1./(epsilon + abs(x-y(i)).^(2-p));
        y1 = sum( w .* x) ./ sum(w);
        ynew(i) = (1-tau)*y(i) + tau*y1;
    end
    plot(transpose([y(:) ynew(:)]), '-', 'color', [s 0 1-s], 'LineWidth', 2);
    y = ynew;
    plot(x, 'k.', 'MarkerSize', 15);
    axis equal; axis tight; axis off;
    drawnow;
    mysaveas(it);
end


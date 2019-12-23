%%
% Use a reweighting l2 scheme for the computation of medians.

n = 300;
p = round(.1*n); % outliers

crandn = @(n)randn(n,1) + 1i*randn(n,1);

x = [crandn(n-p); ...
        crandn(p)*.3 + 30 *( .4+1i*.6 ); ...
        crandn(p)*.3 + 30 *( .6+1i*.6 ); ...
        crandn(p)*.3 + 30 *( .4+1i*.8 )];
y = mean(x);


niter = 50; 
tau = .2; % slow down
epsilon = 1e-9;
for it=1:niter
    clf; hold on;
    plot(x, '.', 'MarkerSize', 15);
    plot(y, '.', 'MarkerSize', 15);
    % axis([-1 1 -1 1]*4);
    axis equal; axis tight;
    drawnow;
    %
    w = 1./(epsilon + abs(x-y));
    y1 = sum( w .* x) ./ sum(w);
    y = (1-tau)*y + tau*y1;
end

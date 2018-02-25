%%
% Test for Metropolis-Hastings for a discrete distribution.

rep = '../results/metropolis-hastings/';
[~,~] = mkdir(rep);

% dimension
n = 30;
% #samples
q = 3000000;

name = 'bimodal';
name = 'triangle';
name = 'unif';

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


niter = 50;

% initial locations
X = ceil(rand(q,1)*n);

% target distribution
mu = rand(n,1);
t = linspace(0,1,n)';
switch name
    case 'triangle'
        mu = [1:n/2, n/2:-1:1]';
    case 'unif'
        mu = abs(t-1/2)<.25;
    case 'bimodal'
        s = .07;
        G = @(c)exp(-(t-c).^2 / (2*s^2));
        mu = .25 + .7*(abs(t-.25)<.1) + 2*(abs(t-.75)<.15);
        mu = .2 + G(.25) + 1.5*G(.75);
end
mu = mu/sum(mu);

% radius for move
r = 3;

disp_mode = 1;

do_save = 1;
save_rate = ceil(niter/50);
save_rate = 1;

for i=1:niter
    nu = hist(X,1:n); nu = nu/sum(nu);
    clf; hold on;
    switch disp_mode
        case 1
            bar(1:n, mu, 1, 'r', 'EdgeColor', 'r');
            bar(1:n, nu, 'b');
            % plot(1:n, mu, 'r', 'LineWidth', 2);
        case 2
            area(1:n, nu);
            plot(1:n, mu, 'r', 'LineWidth', 2);
    end
    axis tight; box on;
    axis([1/2 n+1/2 0 max(mu)*1.03]);
    SetAR(1/2); 
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    if do_save % && mod(i,save_rate)==1
        saveas(gcf, [rep name '-' num2str(i) '.png'], 'png');
    end
    % proposal move
    V = ceil(rand(q,1)*(2*r+1))-r-1;
    X1 = mod(X + V-1,n)+1;
    % accept/reject
    R = mu(X1) ./ mu(X);
    U = rand(q,1)>R;
    X = X1.*(1-U) + X.*U;
end

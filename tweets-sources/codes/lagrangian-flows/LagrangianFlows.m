%%
% Test for Folker-Planck equation discretization using Lagrangian gradient
% flows.

addpath('../toolbox/');
rep = MkResRep();

flowmode = 'stoch';
flowmode = 'partic';


% number of particles
n = 20;
% number of neighbors for entropy estimation
k = 5;

% initial configuration
x = .4 * (2*rand(2,n)-1);
x(1,1:end/2) = x(1,1:end/2) + .5;
x(2,1:end/2) = x(2,1:end/2) + .5;
%
x(1,end/2+1:end) = x(1,end/2+1:end) -.3;
x(2,end/2+1:end) = x(2,end/2+1:end) - .5;

% potential function gradient
gP =  @(x)x;

% weight in front of entropy

kappa = .2; % langevin
m = 20;  % memory for display

kappa = .1; % particle
m = 40;  % memory for display


mu = .01;

tau = .01/2;
tau = .01/3;

niter = 400;
g = [];
r = 0;
X = repmat(x, [1 1 m]);
for i=1:niter
    switch flowmode 
        case 'stoch'
            % SPDE - Lagenvin
            x = x - tau*gP(x) - sqrt(tau)*kappa*randn(2,n);
        case 'partic'
            % NN-computations
            D = distmat(x,x);
            [Di,I] = sort(D, 1);
            I = I(2:k+1,:); Di = Di(2:k+1,:);
            % gradient for entropy  log(|x-xi|^2) -> (x-xi)/|x-xi|^2
            for s=1:2
                xs = x(s,:); xs = xs(:);
                h = repmat(x(s,:), [k 1]) - xs(I); h = h ./ (Di+mu).^2;
                g(s,:) = mean( h );
            end
            % gradient overall
            G = gP(x) - kappa * g;
            % gradient step
            x = x - tau*G;
    end
    % display
    t = (i-1)/(niter-1);
    col = [t 0 1-t];
    %
    X = cat(3, X(:,:,2:end),x);
    clf; hold on;
    PlotBlend(X,col);
    axis equal; 
    axis([-1 1 -1 1]);
    axis off;
    drawnow;
    if mod(i,4)==1
        r = r+1;
        saveas(gcf, [rep flowmode '-' num2str(n) '-' znum2str(r,3) '.png']);
    end
end

% AutoCrop(rep, 'flow-');

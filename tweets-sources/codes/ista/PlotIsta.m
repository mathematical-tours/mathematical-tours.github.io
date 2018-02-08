%%
% Display ista in 2D

rep = '../results/ista/';
[~,~] = mkdir(rep);
% addpath('../toolbox/');


% dictionary, with 2 atoms
D = randn(2,2);

if not(exist('eta'))
    eta = 4;
end
D = [[1;0], [eta;1]];

% thought after solution 
u0 = [1;0];
y = D*u0;

lambda = .3;

a = 1.5; 
q = 501;
t = linspace(-a,a,q);
[Y,X] = meshgrid(t,t);
Z = [X(:),Y(:)]';

mynorm = @(x)norm(x(:));
f = @(Z)1/2*sum(repmat(y, [1 size(Z,2)]) - D*Z).^2 + lambda * sum(abs(Z));

F = reshape(f(Z), [q q]);
% draw
r = 15;
clf; hold on;
imagesc(t,t,F');
% contour(t,t,D',r, 'k--');
colormap(parula(r+1));
axis([-a a -a a]); axis equal; axis off;
% plot the axes
plot([-a a], [0 0], 'k--', 'LineWidth', 2);
plot([0 0], [-a a],  'k--', 'LineWidth', 2);

if not(exist('P')) && 0
    P = [];
    while true
        [x0,x1,button] = ginput(1);
        if button==3
            break;
        end
        P(:,end+1) = [x0;x1];
        plot(x0, x1, 'b.', 'MarkerSize', 25);
    end
end

r = 12;
u = .5/r + (0:r-1)/r;
P = .95*a*[cos(2*pi*u); sin(2*pi*u)];

tau = .05/norm(D)^2;
niter = 1000;
for i=1:size(P,2)
    x = P(:,i);
    X = ista(D,y,x,lambda,tau, niter);
    PlotTrajectory(X(1,:),X(2,:));
    plot(x(1), x(2), 'b.', 'MarkerSize', 25);
    plot(X(1,end), X(2,end), 'r.', 'MarkerSize', 25);
    drawnow;
end
saveas(gcf, [rep 'ista-' num2str( round(100*eta) ) '.png'], 'png');

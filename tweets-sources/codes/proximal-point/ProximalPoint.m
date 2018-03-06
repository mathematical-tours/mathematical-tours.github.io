%%
% test for the proximal point in the quadratic case.

rep = '../results/proximal-point/';
[~,~] = mkdir(rep);
addpath('../toolbox/');


if not(exist('eta'))
    eta = 10;
end

R = @(t)[cos(t), sin(t); -sin(t), cos(t)];
t = .3;
A = R(t)*diag([1,eta])*R(-t);

x0 = [.8;.6];
y = A*x0;

% grid for vizualization
a = 1; 
q = 501;
t = linspace(-a,a,q);
[Y,X] = meshgrid(t,t);
Z = [X(:),Y(:)]';
% function 1/2*<A*z,z> - <z,y>
F = 1/2 * sum( (A*Z).*Z ) - y(1)*Z(1,:) - y(2)*Z(2,:);
F = reshape(F, [q q]);
F = sqrt( F-min(F(:)) );

% draw
r = 20;
clf; hold on;
imagesc(t,t,F');
contour(t,t,F',r, 'k-');
plot(x0(1),x0(2), 'r.', 'MarkerSize', 25);
colormap(parula(r+1));
axis([-a a -a a]); axis equal; axis off;

% proximal point
tau_list = [.01 .1 1 10 100];
for k=1:length(tau_list)
    c = (k-1)/(length(tau_list)-1);
    tau = tau_list(k);
    niter = 300;
    x = [-.5 -.8]';
    for i=1:niter
        x(:,end+1) = (eye(2)+tau*A) \ (x(:,end) + tau*y);
    end
    plot(x(1,:), x(2,:), '.-', 'LineWidth', 2, 'MarkerSize', 25, 'Color', [c 0 1-c]);
end
saveas(gcf, [rep 'proxpt-' num2str(eta) '.png'], 'png');
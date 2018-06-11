%%
% Test for stable marriage

addpath('../toolbox/');
rep = MkResRep();

n = 50;

X = randn(2,n)*.3;

% Y = randn(2,n);
theta = 2*pi*rand(1,n);
r = .8 + .2*rand(1,n);
Y = [cos(theta).*r; sin(theta).*r];

D = distmat(X,Y); % D(xi,yj)

% preference x->y
[D1,A1] = sort(D, 2, 'ascend');
% preference y->x
[D2,A2] = sort(D, 1, 'ascend'); A2 = A2';



S = galeshapley(n, A1, A2); % Men->Women
clf; hold on;
plot([X(1,S); Y(1,:)], [X(2,S); Y(2,:)], 'b-', 'LineWidth', 2);
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 25);
plot(Y(1,:), Y(2,:), 'r.', 'MarkerSize', 25);
axis equal; axis tight;
axis off;
saveas(gcf, [rep 'stable-x-y.eps'], 'epsc');

T = galeshapley(n, A2', A1'); % Women->Men
clf; hold on;
plot([X(1,:); Y(1,T)], [X(2,:); Y(2,T)], 'r-', 'LineWidth', 2);
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 25);
plot(Y(1,:), Y(2,:), 'r.', 'MarkerSize', 25);
axis equal; axis tight;
axis off;
saveas(gcf, [rep 'stable-y-x.eps'], 'epsc');


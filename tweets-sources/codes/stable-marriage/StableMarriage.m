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

I = -ones(n,1);
J = -ones(n,1); % reverse permutation

K = find(I<0);
while not(isempty(K))
    i = K(1); 
    Di = D(i,:); 
    % L = find(not(isinf(Di)));
    [~,j] = min(Di); D(i,j) = +Inf; 
    m = J(j);
    dodraw = 1;
    if m==-1
        I(i) = j; J(j) = i; 
    elseif D(i,j)<D(m,j)
        I(i) = j; J(j) = i;
        I(m) = -1; 
    else
        dodraw = -1;
    end
    % display
    if dodraw     
        f = find(I~=-1);
        clf; hold on;
        plot([X(1,f); Y(1,I(f))], [X(2,f); Y(2,I(f))], 'b-', 'LineWidth', 2);
        plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 25);
        plot(Y(1,:), Y(2,:), 'r.', 'MarkerSize', 25);
        axis equal; axis tight;
        axis off;
        drawnow;
    end
    K = find(I<0);
end

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


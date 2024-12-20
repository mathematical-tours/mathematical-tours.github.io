%  Forward Stagewise linear regression
% Matching pursuits with time-frequency dictionaries
%
%  H. Friedman and W. Stuetzle, "Projection pursuit regression," J, Amer, Statist, Asso., voL 76, pp. 817-823, 1981.
% Friedman and Tukey a projection pursuit algorithm for exploratory data analysis
%
% L K. Jones, "On a conjecture of Huber concerning the convergence of projection pursuit regression," Ann, Statist" voL 15, no. 2, pp, 880-882, 1987

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);





% vectors in the basis
d = 2;
p = 3;
D = randn(d,p);




theta = [pi/3 2*pi/3];
D = [[1;0], [cos(theta(1));sin(theta(1))], [cos(theta(2));sin(theta(2))]];

theta = pi/4;
D = [[1;0], [cos(theta);sin(theta)]];


D = D./sqrt( sum(D.^2,1) );

Col = distinguishable_colors(p+1);
Col(4,:) = []; % remove black

q = 100;
for it=1:q

    theta1 = it/q * 2*pi;
    x0 = [cos(theta1);sin(theta1)];
    
    niter = 70;
    x = zeros(d,1);
    X = x; I = [];
    for jt=1:niter
        E(jt) = norm(x-x0);
        C = D'*(x-x0); % /norm(x-x0);
        [~,i] = max(abs(C));
        x = x - C(i)*D(:,i);
        X(:,end+1) = x;
        I(end+1) = i;
    end
    
     % plot(log10(E));
    
    clf; hold on;
    for s=1:size(X,2)-1
        i = I(s);
        plot(X(1,s:s+1), X(2,s:s+1), '-', 'LineWidth', 2, ...
            'MarkerSize', 15,'color', Col(i,:));
    end
    plot(x(1), x(2), 'k.', 'MarkerSize', 20);
    for k=1:size(D,2)
        plot([0,D(1,k)], [0,D(2,k)], ':', 'LineWidth', 3, 'color', Col(k,:));
    end
    axis equal; 
    axis([-1 1 -1 1]*1.2);
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end
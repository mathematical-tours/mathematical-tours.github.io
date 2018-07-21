%%
% Random graphs
% stochastic block model

%  Paul W. Holland, Kathryn Blackmond Laskey, Samuel Leinhardt. Stochastic blockmodels: first steps. Social
% Networks, 5(2):109?137, 1983.

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% #points
n = 300; 
n = 30; 

% within clusters
p = .2;
% outside clusters
q = .2;

[J,I] = meshgrid(1:n,1:n);
U = rand(n); 
U = U .* (I>J) + U' .* (J>I);
U = U + diag(ones(n,1)+Inf); % remove diagonal

T = @(p,q)[ ones(n/2)*p,ones(n/2)*q; ones(n/2)*q, ones(n/2)*p ];

Q = 3; 
pmax = 10/n; pmin = 10/n;
plist = linspace(pmax,pmin,Q);
r = 2;
qlist = linspace(pmax^(1/r),(.1/n)^(1/r),Q).^r;

Z0 = randn(n,2); 

% Z0(1:end/2,1) = Z0(1:end/2,1)+6;
plist = linspace(30/n,3/n,Q);
qlist = plist;
for i=1:Q
    t = (i-1)/(Q-1);
    %
    p = plist(i); q = qlist(i);
    A = (U<=T(p,q));
    L = - A + diag(sum(A,1));
    %
    [V,D] = eig(L); D = diag(D);
    if 0
    Z = V(:,2:3);  % 2D embedding
    if i>1
        Z = Z * diag(sign(diag(Z1'*Z)));
    end
    Z1 = Z;
    else
        Z = Z0;
    end
    [I,J] = find(sparse(A)); K = find(I>J);
    I = I(K); J = J(K);
    %
    clf; hold on;
    for m=1:length(I);
        plot([Z(I(m),1) Z(J(m),1)], [Z(I(m),2) Z(J(m),2)], '-', 'LineWidth', 1, 'color', [t 0 1-t]);
    end
    plot(Z(:,1), Z(:,2), '.', 'MarkerSize', 25, 'color', 'k');
    axis equal;
    % axis([-1 1 -1 1]*.4);
    axis tight;
    axis off;
    drawnow;
    saveas(gcf, [rep 'graph-' znum2str(i,2) '.png']);
    %

    if 0
    clf;
    bar(1:k, D(2:k+1), 'FaceColor', [t 0 1-t],'EdgeColor', [t 0 1-t]);
    axis tight;
    SetAR(1/2);
    drawnow;
    saveas(gcf, [rep 'eigen-' znum2str(i,2) '.png']);
    end
    % drawnow;
end
% AutoCrop(rep, ['eigen-']);
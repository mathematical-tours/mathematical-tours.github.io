%%
% Clustering via spectral cut.

name = 'annulus';
name = 'twogauss';

addpath('../toolbox/');
rep = MkResRep(name);


n = 400; 

%%
% Generate two clusters

r = 2;
X0 = randn(n,2);


K = 7; % #NN

T0 = rand(n,1)*2*pi;
R0 = randn(n,1);

q = 50;
slist = linspace(0,1,q);

dprev = randn(n,1);

for it=1:q
    
    s = slist(it);
    
    switch name
        case 'twogauss'
            r = 4*(1-s);
            X = X0*.7;
            X(1:n/2,:) = X(1:n/2,:)-r/2;
            X(n/2+1:end,:) = X(n/2+1:end,:)+r/2;
        case 'annulus'
            r = s*.8;
            k0 = round(n*.3);
            R = [R0(1:k0)*(.2+r);  3 + R0(k0+1:end)*(.1+r)]; 
            X = [cos(T0).*R, sin(T0).*R];
            clf; plot(X(:,1), X(:,2), '.', 'MarkerSize', 20);
    end
    % clf; plot(X(:,1), X(:,2), '.');
    
    %%
    % NN graph
    
    %
    D = distmat(X',X');
    D = D + diag(zeros(n,1)+Inf);
    
    [DS,J] = sort(D, 'ascend');
    J = J(1:K,:); J = J(:);
    I = repmat((1:n), [K 1]); I = I(:);
    A = sparse(I,J,ones(K*n,1));
    A = max(A,A');
    
    
    [I,J] = find(sparse(A)); K1 = find(I>J);
    I = I(K1); J = J(K1);
    
    % Laplacian
    L =  diag(sum(A)) - A;
    [V,D] = eig(full(L)); D = diag(D);
    % fiedler vector
    d = V(:,2); d = d-mean(d);
    % re-align according to diagonal direction
    si = sign(sum(d .* dprev));
    if si==0
        si=1; warning('problem');
    end
    d = d * si;
    dprev = d;
    d = rescale(d);
    
    % color
    Col = [d,d*0,1-d];
    clf; hold on;
    for m=1:length(I)
        plot([X(I(m),1) X(J(m),1)], [X(I(m),2) X(J(m),2)], '-', 'LineWidth', 1, 'color', 'k');
    end
    s = ones(n,1)*30; % size
    scatter( X(:,1), X(:,2), s, Col, 'filled' );
    axis equal;
    axis([-1 1 -1 1]*5);
    % axis tight; 
    axis off;
    drawnow;
    
    saveas(gcf, [rep 'nngraph-' znum2str(it,2) '.png']);
    
end

% AutoCrop(rep, 'nngraph-');
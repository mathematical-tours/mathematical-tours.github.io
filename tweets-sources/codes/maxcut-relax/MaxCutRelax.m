%%
% Delorme Poljak ?90
% Goemans Williamson ?94

addpath('../toolbox/');
rep = MkResRep();

name = 'knn';
name = 'regular';
name = 'triangulation';

n = 200; 


X0 = randn(n,2);
r = 0; % separation
P = X0*.7;
P(1:n/2,:) = P(1:n/2,:)-r/2;
P(n/2+1:end,:) = P(n/2+1:end,:)+r/2;



switch name
    case 'knn'
        % NN graph
        K = 10; % #NN
        D = distmat(P',P');
        D = D + diag(zeros(n,1)+Inf);
        %
        [DS,J] = sort(D, 'ascend');
        J = J(1:K,:); J = J(:);
        I = repmat((1:n), [K 1]); I = I(:);   
    case 'triangulation'
        % triangulation
        DT = delaunayTriangulation(P);
        T = DT.ConnectivityList;
        I = T(:);
        J = [T(:,2); T(:,3); T(:,1)];        
    case 'regular' 
        % Uniform grid
        n0 = round(sqrt(n)); n = n0*n0;
        t = linspace(0,1,n0);
        [y,x] = meshgrid(t,t);
        P = [x(:) y(:)];
        %
        U = reshape(1:n, [n0 n0]);
        resh = @(a)a(:);
        I = [resh(U(:,1:end-1)); resh(U(1:end-1,:))];
        J = [resh(U(:,2:end)); resh(U(2:end,:))];    
end

A = sparse(I,J,ones(length(I),1),n,n);
A = max(A,A');


% Laplacian
L =  diag(sum(A)) - A;
L = full(L);

clf; hold on;
gplot(A,P, 'k');
plot(P(:,1), P(:,2), 'b.', 'MarkerSize', 25);
axis equal; axis tight; axis off;

cvx_solver sdpt3 
%  SeDuMi %
cvx_begin sdp % quiet
cvx_precision high;
variable X(n,n) symmetric;
minimize sum(X(:).*A(:)); % -sum(X(:).*L(:));
subject to
    diag(X)==ones(n,1);
    X >= 0;
cvx_end

% X = U*diag(S)*V' = Z'*Z`

[U,S] = eig(X); S = max(0,diag(S));
U = real(U);
Z = diag(sqrt(S)) * U'; % each column is a vector

r = randn(n,1);
s = Z'*r;

clf; hold on;
gplot(A,P, 'k');
plot(P(s>0,1), P(s>0,2), 'b.', 'MarkerSize', 25);
plot(P(s<=0,1), P(s<=0,2), 'r.', 'MarkerSize', 25);
axis equal; axis tight; axis off;
saveas(gcf, [rep 'maxcut-' name '.eps'], 'epsc');

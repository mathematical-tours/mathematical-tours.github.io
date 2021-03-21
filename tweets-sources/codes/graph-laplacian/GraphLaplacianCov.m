%%
% Display covariance induced by graph laplacian inverse

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 300;
z = rand(n,1)+1i*rand(n,1);
z = 2*z-1-1i;

% uniform graph
if 0
n0 = 20;
n = n0*n0;
x = linspace(-1,1,n0)'+randn(n0,1)*1e-6;
[Y,X] = meshgrid(x,x);
z = X(:)+1i*Y(:);
end

% k-nn graph
k = 4;
D = abs(z-transpose(z));
[D1,I] = sort(D,2, 'ascend'); I = I(:,2:k+1);
J = repmat((1:n)',[1 k]);
A = sparse(I(:),J(:),ones(n*k,1), n, n);
A = max(A,A');

% delaunay graph
DT = delaunayTriangulation([real(z), imag(z)]);
I = [DT.ConnectivityList(:,1);DT.ConnectivityList(:,2);DT.ConnectivityList(:,3)];
J = [DT.ConnectivityList(:,2);DT.ConnectivityList(:,3);DT.ConnectivityList(:,1)];
A = sparse(I(:),J(:),ones(length(I(:)),1));
  


L = full(diag(sum(A))-A);

lambda = .0001;
lambda = 1e-4;
C = inv(L+eye(n)*lambda);

t = 100;
% C = expm(-t*L);


q = 75;
for it=1:q
    s = (it-1)/(q-1);  
    z0 = exp(2i*pi*s)*.7;
    [~,i0] = min(abs(z-z0)); 
    f = rescale(C(i0,:));    
    clf; disp_knn_graph(I,J,z,f);
    scatter(real(z(i0)),imag(z(i0)), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth',2);
    axis off;
    drawnow;
    mysaveas(it);
end
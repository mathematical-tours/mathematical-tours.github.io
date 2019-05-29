%%
% Compute the page-rank eignvector for a graph.

addpath('../toolbox/');
rep = MkResRep();


n = 200; 
% euclidean NN graph

x = randn(n,1) + 1i*randn(n,1);
D = abs(x*ones(1,n)-ones(n,1)*transpose(x));
if 0
    d = sort(abs(D(:)));
    epsilon = d(10*n);
    A = D<epsilon;
    A = A-diag(diag(A));
    [i,j] = find(sparse(A));
    [i,j] = deal(i(i<j), j(i<j));
else
    K = 5; % #NN
    [D1,I] = sort(D, 'ascend');
    j = I(1:K,:)'; j = j(:);
    i = repmat((1:n)',[1 K]); i = i(:);
    A = sparse(i,j,ones(K*n,1));
end



% i = (1:n-1)'; j = (2:n)';
% A = sparse([i;j],[j;i],ones(2*n-2,1));

clf; hold on;
plot(transpose([x(i) x(j)]), 'k', 'LineWidth', 2);
plot(x, 'k.', 'MarkerSize', 25);

% iteration of message passing
u = rand(n,1);
u = ones(n,1)/n;
niter = 50;
M = max(.00001,sum(A,2));
tau = .1;
d = .95;
for it=1:niter
    % u = u/sum(u);
    % display
    clf; hold on;
    plot(transpose([x(i) x(j)]), 'k', 'LineWidth', 2);
    % 2D plot
    % s = .01 + rescale(u)*100; % size
    s = u*n*120;
    sc = rescale(s);
    sc = min(u*n*.7,1);
    Col = sc*[1 0 0] + (1-sc)*[0 0 1];
    scatter( real(x), imag(x), s, Col, 'filled' );
    axis off;
    axis equal; axis tight;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    % axis equal; axis tight; box on;
    u1 = (1-d)/n + d * A'*(u ./ M(:));
    u = (1-tau)*u + tau*u1;
end


%  AutoCrop(rep, 'anim')
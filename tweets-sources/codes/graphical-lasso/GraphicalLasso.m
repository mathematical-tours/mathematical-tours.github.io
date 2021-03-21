
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


%%
% Graphical lasso.
% From 
%   Friedman et al, ?Sparse inverse covariance estimation with the graphical lasso?, 
%   Biostatistics 9, pp 432, 2008
%       --> improve block coord (lasso type)

% REF:  
%  - d'aspremont 1st order https://arxiv.org/pdf/math/0609812.pdf
%       --> first to propose block coord algo
%  - Appli irmf : https://hal.inria.fr/file/index/docid/530853/filename/paper.pdf
%  - spice algo https://www.stat.berkeley.edu/~bickel/Rothman%20et%20al%202008.pdf

% Meinshausen, N. & Buhlmann, ¨ P. (2006), ?High dimensional graphs and variable selection with the lasso?, Annals of Statistics 34, 1436?1462.
% --> first proposal of sparse covariance

% https://people.eecs.berkeley.edu/~elghaoui/Pubs/CvxTechCovSel_ICML.pdf
% + data genomique

dotp = @(x,y)sum(x(:).*y(:));

p = 30; 

% sparse precision matrix

% random graph
mysym = @(x)(x+x')/2;
U = rand(p)<.05; U = max(U,U');
V = mysym(1+rand(p));
V = V.*sign(mysym(randn(p)));
P0 = U .* V; 

% generate sparse k-nn euclidean graph
k = 3; 
x = sort(randn(p,1)) + 1i*randn(p,1);
[D,I] = sort(abs(x-transpose(x)),2); I = I(:,2:k+1);
J = repmat((1:p)', [1 k]);
A = sparse(I,J,ones(p*k,1), p,p); A = max(A,A');

clf; 
plot(graph(A),'XData',real(x),'YData', imag(x), 'LineWidth', 2, 'MarkerSize', 10);
axis equal;

M = 1+rand(p); 
P0 = full(A) .* max(M,M');
P0 = P0 - min(eig(P0))*eye(p);
P0 = P0/max(eig(P0)) + .1*eye(p);

C0 = inv(P0);

% remove diagonal
rem_diag = @(C)C-diag(diag(C));

clf;
subplot(1,2,1);
imagesc(rem_diag(P0)); axis image;
subplot(1,2,2);
imagesc(rem_diag(C0)); axis image;

% sample from the Gaussian 
n = 1000;
X = ( sqrtm(C0)*randn(p,n) )';
% empirical covariance
C_hat = X'*X/n;

clf;
subplot(1,2,1);
imagesc(rem_diag(C0)); axis image;
subplot(1,2,2);
imagesc(rem_diag(C_hat)); axis image;

clf
plot([sort(eig(C_hat)), sort(eig(C0))]);

% solve
%     min_P <P,C> - logdet(P) + lambda*|P|_1


% solve dual problem
%  p + max_{ |C-C_hat|_inf<=lambda } logdet(C)

lambda = .3;
[C,E] = glasso(C_hat,lambda);
P = inv(C);

clf;
subplot(1,2,1);
imagesc(rem_diag(P0));
subplot(1,2,2);
imagesc(rem_diag(P));

clf;
plot(log10(E(1:end/2)-min(E)));

lambda_list = [.01 .1 .2 .5];
clf; 
for it=1:length(lambda_list)
    lambda = lambda_list(it);
    [C,E] = glasso(C_hat,lambda);
    P = inv(C);
    A1 = sparse(abs(P)>1e-2);
    subplot(2,2,it); hold on;
%    plot(graph(A, 'omitselfloops'),'k-','XData',real(x),'YData', imag(x), 'LineWidth', 1, 'MarkerSize', 3,'NodeLabel',{});
    plot(graph(A1, 'omitselfloops'),'r','XData',real(x),'YData', imag(x), 'LineWidth', 2, 'MarkerSize', 3,'NodeLabel',{});
    axis equal; axis off;
end

q = 50; 
lambda_list = linspace(.05,1,q); 
for it=1:length(lambda_list)
    lambda = lambda_list(it);
    [C,E] = glasso(C_hat,lambda);
    P = inv(C);
    A1 = sparse(abs(P)>1e-2);
	clf;  hold on;
    % plot(graph(A, 'omitselfloops'),'k-','XData',real(x),'YData', imag(x), 'LineWidth', 1, 'MarkerSize', 3,'NodeLabel',{});
    plot(graph(A1, 'omitselfloops'),'r','XData',real(x),'YData', imag(x), 'LineWidth', 2, 'MarkerSize', 3,'NodeLabel',{});
    axis equal; axis off;
    mysaveas('graph', it);
    
    clf;
    imagesc(rem_diag(P));
    axis image; axis off;
    caxis( [0,max(max(rem_diag(P0)))] );
    mysaveas('precis', it);
    
    clf;
    imagesc(rem_diag(C));
    axis image; axis off;
    caxis( [0,max(max(rem_diag(C0)))] );
    drawnow;
    mysaveas('cov', it);
end


return;

% compare with another solver
[W1,C1] = glasso_solver(C_hat, lambda);

plot(C1(:)-C_hat(:), '.')

clf; hold on;
plot(sort(eig(C)), 'b');
plot(sort(eig(C1)), 'r--');

return;

lambda = .2;
[P,W] = glasso_solver(hatC, lambda);

clf; 
plot([sort(eig(P)), sort(eig(P0))]);
axis([1 p 0 max(eig(P0))*1.1]);

s = max(max(rem_diag(P0)))
clf;
subplot(1,2,1);
imagesc(rem_diag(P)); axis image;
caxis([-.2*s,s]);
subplot(1,2,2);
imagesc(rem_diag(P0)); axis image;
caxis([-.2*s,s]);



return;



thresh = @(x,t)max(abs(x)-t,0).*sign(x);
threshD = @(P,t)thresh(P-diag(diag(P)),t)+diag(diag(P));

lambda = .0;

niter = 1000;
tau = .08/10;
tau_max = .1; ntau = 40;
Elist = [];
P = eye(p);
for it=1:niter
    [Elist(end+1),gradE] = glasso_energy(P,C,lambda);
    % linesearch
    if 0
        e = [];
        for j=1:ntau
            tau = tau_max*(j-1)/(ntau-1);
            Q = threshD(P-tau*g, tau*lambda);
            e(j) = glasso_energy(Q,C,lambda);
        end      
        [~,j] = min(e); tau = tau_max*(j-1)/(ntau-1);
        if j==1
            warning('not progressing');
            tau_max = tau_max/5;
            % return;
        end
    end
    %
    P = threshD(P-tau*g, tau*lambda);
    clf; 
    plot([sort(eig(P)), sort(eig(P0))]);
    axis([1 p 0 max(eig(P0))*1.1]);
    drawnow;
end

clf;
subplot(1,2,1);
imagesc(rem_diag(P)); axis image;
subplot(1,2,2);
imagesc(rem_diag(P0)); axis image;
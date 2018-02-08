%%
% Bench for Lostitic classification.

rep = '../results/logistic/';
[~,~] = mkdir(rep);
addpath('../toolbox/');


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);
dotp = @(u,v)sum(u(:).*v(:));
AddBias = @(X)[X ones(size(X,1),1)];

% click selection
if not(exist('Z'))
    it = 0;
    clf; hold on;
    Z = [];
    while true
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if button==3
            break;
        end
        Z(:,end+1) = [a;b];
    end
end

k = size(Z,2); % #classes
% width of classes
if not(exist('sigma'))
    sigma = .05;
end

n0 = 500; % number of sample per class
n = n0*k;
p = 2; % dimensionality
X = []; y = [];
randn('state', 123);
sigma0 = .05;
for i=1:k
    U = randn(n0,2)*sigma0;
    z = Z(:,i); z = (z-.5)*sigma+.5;
    U(:,1)=U(:,1)+z(1); U(:,2)=U(:,2)+z(2);
    X = [X;U];
    y = [y; i*ones(n0,1)];
end
% proba matrix
CL = 1:k;
D = double( repmat(CL(:)', [n,1]) == repmat(y, [1,k]) );

% LSE and SM
X1 = AddBias(X);
LSE = @(S)log( sum(exp(S), 2) );
max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);
LSE = @(S)LSE( S-max2(S) ) + max(S,[],2);
SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);
SM = @(S)SM(S-max2(S));
% energy
E = @(W)1/n*( sum(LSE(X1*W)) - dotp(X1*W,D)  );
nablaE = @(W)1/n * X1'* ( SM(X1*W) -  D );

% gradient descent
W = zeros(p+1,k);
Elist = [];
tau = 2;
niter = 8000;
for i=1:niter
    W = W - tau * nablaE(W);
    Elist(i) = E(W);
end
clf; plot(Elist);
if Elist(end)>Elist(1)
    warning('Convergence problem');
end

% grid to evaluate the class density
M = .03;
q = 201;
t = linspace(-M,1+M,q);
[B,A] = meshgrid(t,t);
G = zeros(q*q,p);
G(:,1:2) = [A(:), B(:)];
%
Theta = SM(AddBias(G)*W);
Theta = reshape(Theta, [q q k]);

% render in color
col = distinguishable_colors(k+1)'; col = col(:,[1:3 5:end]);
R = zeros(q,q,3);
for i=1:k
    for a=1:3
        R(:,:,a) = R(:,:,a) + Theta(:,:,i) .* col(a,i);
    end
end
% display
clf; hold on;
imagesc(t,t,permute(R, [2 1 3]));
for i=1:k
    I = find(y==i);
    plot(X(I,1),X(I,2), '.', 'color', col(:,i)*.6, 'MarkerSize', 21);
end
axis equal; axis([-M,1+M,-M,1+M]); axis off;
saveas(gcf, [rep 'logistic-k' num2str(k) '-' num2str(round(100*sigma)) '-regul.png'], 'png');

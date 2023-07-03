%%
% Test of adaboost for axis-aligned classifiers.

rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);
mysaveas = @(it)0;

%%
% generate data sets

n = 1000; p = 2;

% cylinder
t = 2*pi*rand(n/2,1);
R = .95;
r = R*(1 + .08*randn(n/2,1)); % radius
X1 = [cos(t).*r, sin(t).*r];
X = [.25*randn(n/2,2); X1];

% two blobs
y = [ones(n/2,1);-ones(n/2,1)];

%%
% Display 

clf; hold on;
plot(X(y>0,1), X(y>0,2), '.r', 'MarkerSize', 20);
plot(X(y<0,1), X(y<0,2), '.b', 'MarkerSize', 20);
axis equal; axis([-1 1 -1 1]);
set(gca, 'Xtick', [], 'Ytick', []); box on;

%%
% Define the weak learners as axis-aligned learning

r = 100;  % number of weak learners


% random orientation
theta = rand(1,r)*2*pi;
u = [cos(theta); sin(theta)];
% random intercept
t = 1.2 * (2*rand(1,r)-1);
% decision
a = 1/n * sum( sign(  (X(:,1)*u(1,:) + X(:,2)*u(2,:)-t) ) == y );
sigma = sign(a-1/2 + 1e-10);
% error, should be >1/2 for weak-learning property
H = sign(  (X(:,1)*u(1,:) + X(:,2)*u(2,:)-t).*sigma  );
E = ( H ~= y );
e = 1/n * sum( E );

%%
% Grid for display

m = 300;
s = linspace(-1,1,m);
[V,U] = meshgrid(s,s);
Z = [U(:),V(:)];

%%
% Display a few weak classifiers ...

for k=1:12
    y1 = sign( (Z(:,1)*u(1,k) + Z(:,2)*u(2,k)-t(k)).*sigma(k)  );
    y1 = reshape(y1,[m,m]);
    
    clf; hold on;
    imagesc(s,s,-y1'); colormap gray(256); caxis([-1 1])
    sz = ones(n,1)*25;
    scatter( X(y>0,1), X(y>0,2),sz(y>0), [1 0 0], 'filled' );
    scatter( X(y<0,1), X(y<0,2),sz(y<0), [0 0 1], 'filled' );
    axis equal; axis([-1 1 -1 1]);
    % axis off;
    set(gca, 'Xtick', [], 'Ytick', []); box on;
    drawnow;
    % saveas(gcf, [rep 'weak-' znum2str(k,2) '.png']);
end


%%
% update

% weights on the model
w = zeros(1,r);
% probability on sample
D = ones(n,1)/n;

niter = 2000;
q = 100;
idisp = 1;
ndisp = round( 1+(niter-1)*linspace(0,1,q).^4 );
ndisp = unique(ndisp);

Err = []; % addboost minimize a convex exponential upper bound on the empirical risk
for it=1:niter
    % evaluate the error of the classifier
    Err(it) = 1/n * sum( exp(-y .* sum(w .* H,2)) );
    % minimum error learner
    [e,i] = min( sum( D .* E  ) );
    % gradient step
    a = 1/2 * log( (1-e)/e );
    w(i) = w(i) + a;
    % proba update update
    D = D .* exp(-a*y.*H(:,i));
    D = D / ( 2*sqrt(e*(1-e)) );
    %D = D/sum(D);
    % evaluate on a grid
    % Hw =  sign((Z(:,d)-t).*sigma); % weak learner eval
    Hw = sign( (Z(:,1)*u(1,:) + Z(:,2)*u(2,:)-t).*sigma  );
    y1 = ( sum(w .* Hw,2) ); % bagging aggregation
    y1 = reshape(y1,[m,m]);
    % display
    if ndisp(idisp)==it
        clf; hold on;
        imagesc(s,s,-y1'); colormap gray(256); caxis([-1 1])
        sz = rescale(D,.3,1)*50;
        scatter( X(y>0,1), X(y>0,2),sz(y>0), [1 0 0], 'filled' );
        scatter( X(y<0,1), X(y<0,2),sz(y<0), [0 0 1], 'filled' );        
        axis equal; axis([-1 1 -1 1]);
        % axis off;
        set(gca, 'Xtick', [], 'Ytick', []); box on;
        drawnow;
        mysaveas(idisp);
        idisp = idisp+1;
    end
end
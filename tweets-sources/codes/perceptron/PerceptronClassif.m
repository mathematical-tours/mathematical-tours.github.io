%%
% Test for classification in 2 classes in 2-D using a perceptron with a single
% hidden layer and a logistic loss.


name = 'twomoons';
name = 'handmade';

addpath('../toolbox/');
rep = MkResRep('classif');



% sigmoid
phi = @(x)1./(1+exp(-x));
phiD = @(x)exp(-x)./(1+exp(-x)).^2;

% soft max
max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);
SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);
SM = @(S)SM(S-max2(S));

% Solvers
% min_{a,b,c} E(A,c) = sum_i |y_i-f(x_i,a,b,c)|^2
% f(x,a,b,c) = sum_i c_i phi(<a_i,x>-b_i)  = phi(XA-b)*c

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

t = linspace(-5,5,1000);
clf; hold on;
plot(phiD(t));
plot((phi(t+1e-5) - phi(t))/(1e-5), 'k--');


n0 = 200; % number of points per class
p = 5;   % number of inner neurons
k = 2;   % number of classes
n = n0*k; % Total number of points

col = {'r' 'g' 'b'};

switch name
    case 'twomoons'
        Bmin = [-1.1 -1.1]; Bmax = [1.1 2.1];
        eta = .13;
        theta = [rescale(rand(n0,1),-eta,1+eta)*pi; pi + rescale(rand(n0,1),-eta,1+eta)*pi]; % angles
        r = 1 + randn(n0*2,1)*.06;
        x = [ r.*cos(theta), r.*sin(theta) ];
        x(end/2+1:end,1) = x(end/2+1:end,1)  + 1;
        y = [ones(n0,1);-ones(n0,1)];
    case 'handmade'
        Bmin = [0 0]; Bmax = [1 1];
        if not(exist('y'))
            [x,y] = PickDistributions();
        end
end

n = size(x,1);

X = [x, ones(n,1)]; % add bias
d = size(X,2); % dimension of the data (2D+offstet)
cl = [-1 1];


clf;  hold on;
for j=1:k
    I = find(y==cl(j));
    plot(x(I,1), x(I,2), '.', 'color', col{j}, 'MarkerSize', 20);
end
axis equal; axis off;


Phi  = @(A)phi(X*A);
PhiD = @(A)phiD(X*A);

% L2 loss
E = @(A,c)1/(2*n)*norm(Phi(A)*c-y)^2;
% grad of loss
R = @(A,c)Phi(A)*c-y;

% L2 loss
Loss = @(s,y)1/2*norm(s-y)^2;
NablaLoss = @(s,y)s-y;

% Logistic loss
Loss = @(s,y)sum( log( 1 + exp(-s.*y) ) );
theta = @(v)1 ./ (1+exp(-v));
NablaLoss = @(s,y)- y.* theta(-s.*y);

% grad of loss
E = @(A,c)1/n * Loss(Phi(A)*c,y);
R = @(A,c)NablaLoss(Phi(A)*c,y);
%
GEA = @(A,c)1/n * X' * ( PhiD(A) .* ( R(A,c)*c' ) );
GEc = @(A,c)1/n * Phi(A)' * R(A,c);




% check finite difference
rho = 1e-6;
nrep = 10;
dotp = @(u,v)sum(u(:).*v(:));
A = randn(d,p); c = randn(p,1);
for i=1:nrep
    %
    delta = randn(d,p);
    f = @(A)E(A,c); g = GEA(A,c);
    h = ( f(A+delta*rho)-f(A) )/rho; 
    Err_A(i) = ( dotp(delta,g) - h ) / h;
    %
    delta = randn(p,1);
    f = @(c)E(A,c); g = GEc(A,c);
    h = ( f(c+delta*rho)-f(c) )/rho; 
    Err_c(i) = ( dotp(delta,g) - h ) / h;
end
clf; 
subplot(2,1,1); plot(Err_A);
subplot(2,1,2); plot(Err_c);


% initialization
A0 = 10*randn(d,p);
c0 = randn(p,1);

% c0 = Phi(A)'*y;
c0 = .5*c0 * norm(y)/norm(Phi(A)*c0);

% intercept for BFGS
a = A0; c = c0;
flat = @(x)x(:);
Aa=@(r)reshape(r(1:d*p),[d p]);  C=@(r)r(d*p+1:end);
GradE = @(A,c)[flat(GEA(A,c));GEc(A,c)];
F = @(r)deal(   E(Aa(r),C(r)), ...
                GradE(Aa(r),C(r))  );
% BFGS
r0 = [A0(:);c0];
options.niter = 400;
options.report = @(f,v)struct('val',v,'rlist',f);
[r, R, info] = perform_bfgs(F, r0, options);
rlist = s2v(R, 'rlist');
Elist = s2v(R, 'val'); 
clf; plot(Elist, 'LineWidth', 2);


% display list
N0 = 100; N = N0*N0;
u = linspace(Bmin(1),Bmax(1),N0);
v = linspace(Bmin(2),Bmax(2),N0);
[V,U] = meshgrid(v,u);
Z = [V(:), U(:), ones(N,1)]; % add bias


lvl = 20; % #levellines
vmax = 1.2; % max value for display


% rendering
niter = min(length(Elist),200);
q = 50;
disp_list = round(linspace(1,niter,q));
it = 1;
for i=1:niter
    if i==disp_list(it)
        r = rlist(:,i); A = Aa(r); c = C(r);
        t = (i-1)/(niter-1);  
        % evaluate on the grid
        M = reshape( theta( phi(Z*A)*c ), [N0 N0]);
        % display
        clf; hold on;
        imagesc(v,u,M);
        contour(v,u,M,linspace(-vmax,vmax,lvl), 'k');
        colormap(parula(lvl-1));
        caxis([-vmax vmax]);
        axis image; axis off;
        %
        for j=1:k
            I = find(y==cl(j));
            plot(x(I,1), x(I,2), '.', 'color', col{j}, 'MarkerSize', 20);
        end
        axis equal; axis off;
        %
        drawnow;
        saveas(gcf, [rep name '-sol-' num2str(p) '-' znum2str(it,3) '.png'], 'png');
        it = it+1;
    end
end

% Rendering of energy
it = 1;
for i=1:niter
    % display
    if i==disp_list(it)
        clf; hold on;
        plot(Elist(1:niter), 'k', 'LineWidth', 2);
        plot(i, Elist(i), '.', 'Color', [t,0,1-t], 'MarkerSize', 50);
        axis([0 niter 0 .7]);
        box on; set(gca, 'FontSize', 15);
        SetAR(1/2);
        drawnow;
        saveas(gcf, [rep name '-ener-' num2str(p) '-' znum2str(it,3) '.png'], 'png');
        it = it+1;        
    end
end

% AutoCrop(rep, [name '-sol-' num2str(p) '-']); 
% AutoCrop(rep, [name '-ener-' num2str(p) '-']); 



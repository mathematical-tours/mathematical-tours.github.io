%%
% Test for approximation of 1D function using a perceptron with a single
% hidden layer.

addpath('../toolbox/');
rep = MkResRep('2d');

% sigmoid
phi = @(x)1./(1+exp(-x));
phiD = @(x)exp(-x)./(1+exp(-x)).^2;

% Solvers
% min_{a,b,c} E(A,c) = sum_i |y_i-f(x_i,a,b,c)|^2
% f(x,a,b,c) = sum_i c_i phi(<a_i,x>-b_i)  = phi(XA-b)*c


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

t = linspace(-5,5,1000);
clf; hold on;
plot(phiD(t));
plot((phi(t+1e-5) - phi(t))/(1e-5), 'k--');

n0 = 80; 
n = n0*n0; % #samples
p = 100;  % #hidden layer

xmax = 2.5; 
u = linspace(-xmax,xmax,n0)';
[V,U] = meshgrid(u,u);
X = [U(:), V(:), ones(n,1)]; % add bias
d = size(X,2); % dimension of the data (2D+offstet)


name = 'sine';
name = 'peaks';
% target function 
switch name
    case 'sine'
        y = ( sin(3*pi*abs(U).^1.5) + abs(U) ) ...
        .*  ( sin(3*pi*abs(V).^1.5) + abs(U) );
    case 'peaks'
        y = peaks(U,V); 
end
y = rescale(y(:), -1,1);

Phi  = @(A)phi(X*A);
PhiD = @(A)phiD(X*A);

resh = @(y)reshape(y,[n0 n0]);

% loss
E = @(A,c)1/(2*n)*norm(Phi(A)*c-y)^2;
% grad
R = @(A,c)Phi(A)*c-y;
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

% render target function
lvl = 20; % #levellines
vmax = 1.25; % max value for display
clf; hold on;
imagesc(u,u,resh(y));
contour(u,u,resh(y),linspace(-vmax,vmax,lvl), 'k');
colormap(parula(lvl-1));
caxis([-vmax vmax]);
axis image; axis off;
drawnow;
saveas(gcf, [rep name '-input.png'], 'png');


% rendering
niter = length(Elist);
ndisp = ceil(niter/300);
for i=1:niter
    if ndisp==1 || mod(i,ndisp)==1
        r = rlist(:,i); A = Aa(r); c = C(r);
        M = resh(Phi(A)*c);
        t = (i-1)/(niter-1);      
        % display
        clf; hold on;
        imagesc(u,u,M);
        contour(u,u,M,linspace(-vmax,vmax,lvl), 'k');
        colormap(parula(lvl-1));
        caxis([-vmax vmax]);
        axis image; axis off;
        drawnow;
        saveas(gcf, [rep name '-sol-' num2str(p) '-' znum2str((i-1)/ndisp+1,3) '.png'], 'png');
    end
end

% Rendering of energy
for i=1:niter
    % display
    if ndisp==1 || mod(i,ndisp)==1
        clf; hold on;
        plot(Elist, 'k', 'LineWidth', 2);
        plot(i, Elist(i), '.', 'Color', [t,0,1-t], 'MarkerSize', 50);
        axis([0 niter 0 .05]);
        box on; set(gca, 'FontSize', 15);
        SetAR(1/2);
        drawnow;
        saveas(gcf, [rep name '-ener-' num2str(p) '-' znum2str((i-1)/ndisp+1,3) '.png'], 'png');
    end
end

% AutoCrop(rep, [name '-sol-' num2str(p) '-']); 
% AutoCrop(rep, [name '-ener-' num2str(p) '-']); 



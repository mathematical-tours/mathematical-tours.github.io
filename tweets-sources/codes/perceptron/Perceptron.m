%%
% Test for approximation of 1D function using a perceptron with a single
% hidden layer.

addpath('../toolbox/');
rep = MkResRep();

% sigmoid
phi = @(x)1./(1+exp(-x));
phiD = @(x)exp(-x)./(1+exp(-x)).^2;

% Solvers
% min_{a,b,c} E(a,b,c) = sum_i |y_i-f(x_i,a,b,c)|^2
% f(x,a,b,c) = sum_i c_i phi(a_i*x-b_i) 


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

t = linspace(-5,5,1000);
clf; hold on;
plot(phiD(t));
plot((phi(t+1e-5) - phi(t))/(1e-5), 'k--');

n = 256; % #samples
p = 6;  % #hidden layer
x = linspace(-1,1,n)';

name = 'bumps';
name = 'sine';
% target function 
switch name
    case 'sine'
        y = sin(6*pi*abs(x).^1.5) + abs(x)*2;
    case 'bump'
        s = .1;
        y = exp(-x.^2/(2*s^2));
    case 'bumps'
        s = .1; c = [-.5 .5];
        y = exp(-x.^2/(2*s^2));
        y = exp(-(x-c(1)).^2/(2*s^2)) - 1.2*exp(-(x-c(2)).^2/(2*s^2)) ;
end
Phi  = @(a,b)phi(x*a'+ones(n,1)*b');
PhiD = @(a,b)phiD(x*a'+ones(n,1)*b');

% loss
E = @(a,b,c)1/(2*n)*norm(Phi(a,b)*c-y)^2;
% grad
R = @(a,b,c)Phi(a,b)*c-y;
GEa = @(a,b,c)1/n * ( PhiD(a,b)' .* ( c * R(a,b,c)' ) ) * x;
GEb = @(a,b,c)1/n * ( PhiD(a,b)' .* ( c * R(a,b,c)' ) ) * ones(n,1);
GEc = @(a,b,c)1/n * Phi(a,b)' * R(a,b,c);

% check finite difference
rho = 1e-6;
nrep = 10;
dotp = @(u,v)sum(u(:).*v(:));
a = randn(p,1); b = randn(p,1); c = randn(p,1);
for i=1:nrep
    delta = randn(p,1);
    %
    f = @(a)E(a,b,c); g = GEa(a,b,c);
    d = ( f(a+delta*rho)-f(a) )/rho; 
    Err_a(i) = ( dotp(delta,g) - d ) / d;
    %
    f = @(b)E(a,b,c); g = GEb(a,b,c);
    d = ( f(b+delta*rho)-f(b) )/rho; 
    Err_b(i) = ( dotp(delta,g) - d ) / d;
    %
    f = @(c)E(a,b,c); g = GEc(a,b,c);
    d = ( f(c+delta*rho)-f(c) )/rho; 
    Err_c(i) = ( dotp(delta,g) - d ) / d;
end
clf; 
subplot(3,1,1); plot(Err_a);
subplot(3,1,2); plot(Err_b);
subplot(3,1,3); plot(Err_c);


% initialization
a0 = 10*ones(p,1)+randn(p,1)*3;
a0 = 10*randn(p,1);
b0 = linspace(-1.2,1.2,p)'; 
b0 = randn(p,1)*2;
c0 = randn(p,1);

% c0 = Phi(a,b)'*y;
c0 = .5*c0 * norm(y)/norm(Phi(a,b)*c0);

% intercept for BFGS
a = a0; b = b0; c = c0;
A=@(r)r(1:p);  B=@(r)r(p+1:2*p);  C=@(r)r(2*p+1:3*p);
GradE = @(a,b,c)[GEa(a,b,c);GEb(a,b,c);GEc(a,b,c)];
F = @(r)deal(   E(A(r),B(r),C(r)), ...
                GradE(A(r),B(r),C(r))  );
% BFGS
r0 = [a0;b0;c0];
options.niter = 400;
options.report = @(f,v)struct('val',v,'rlist',f);
[r, R, info] = perform_bfgs(F, r0, options);
rlist = s2v(R, 'rlist');
Elist = s2v(R, 'val'); 
clf; plot(Elist, 'LineWidth', 2);

% rendering
niter = length(Elist);
ndisp = ceil(niter/300);
for i=1:niter
    % display
    if ndisp==1 || mod(i,ndisp)==1
        r = rlist(:,i); a = A(r); b = B(r); c = C(r);
        t = (i-1)/(niter-1);
        clf; hold on;
        plot(x, y, 'k', 'LineWidth', 2);
        plot(x, Phi(a,b)*c, 'LineWidth', 2, 'Color', [t,0,1-t]);
        axis([-1,1,min(y)-.5,max(y)+.5]); 
        box on; set(gca, 'XTick', [], 'YTick', []);
        SetAR(2/3);
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
        axis([0 niter 0 .5]);
        box on; set(gca, 'FontSize', 15);
        SetAR(1/2);
        drawnow;
        saveas(gcf, [rep name '-ener-' num2str(p) '-' znum2str((i-1)/ndisp+1,3) '.png'], 'png');
    end
end

% AutoCrop(rep, [name '-sol-' num2str(p) '-']); 
% AutoCrop(rep, [name '-ener-' num2str(p) '-']); 



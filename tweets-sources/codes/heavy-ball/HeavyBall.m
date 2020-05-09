%%
% Test for first/second order ODE for function minimization.

n = 256; % for display

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
Z = X + 1i*Y;

c = 1.1*[-.5-.2i, +.5-.2i, .1+.4i];
w = [1 .8 .9 .6];
s = .2; 
phi = @(z,c)1-exp(-abs(z-c).^2/(2*s^2));
f = @(z) ...
    w(1)*phi(z,c(1)) + ...
    w(2)*phi(z,c(2)) + ...
    w(3)*phi(z,c(3)) + w(4)*abs(z).^2;
nablaphi = @(z,c)(z-c)/(s^2) .* exp(-abs(z-c).^2/(2*s^2));
nablaf = @(z) ...
    w(1)*nablaphi(z,c(1)) + ...
    w(2)*nablaphi(z,c(2)) + ...
    w(3)*nablaphi(z,c(3)) + 2*w(4)*z;


K = 5;  % #particles
r = 12; % #levellines
m = linspace(0,1,r-1)';
CM = m*[1 1 1] + (1-m)*[.3 .3 .3];
col = distinguishable_colors(K);

% seed points
x0 = .98*exp( 2i*pi * (0:K-1)'/K );





% gradient descent
alpha = Inf; 
tau = .01/4; 
tmax = .9;




% momentum/heavyball
nest_mode = 0; 
alpha = .2;  
tau = .05; 
tmax = 10;


% momentum/heavyball
nest_mode = 0; 
alpha = 1;  
tau = .05; 
tmax = 10;


% nesterov 
nest_mode = 1; 
alpha = 1; 
tau = .05; 
tmax = 10;


% momentum/heavyball
nest_mode = 0; 
alpha = 0;  
tau = .05; 
tmax = 10;

niter = tmax/tau;

x = x0; 
v = x0*0; % speed
for j=1:niter
    switch alpha
        case Inf
            % pure gradient descent: x' + nablaf(x)=0
            x(:,end+1) = x(:,end) - tau*nablaf( x(:,end) );
        otherwise
            if nest_mode==1
                a = 3/( j*tau );
            else
                a = alpha;
            end
            % momentum: x'' + alpha*x' + nablaf(x)=0
            % x'=v,  v' = -alpha*v - nablaf(x)=0  ---> leapfrog integrator
            v = v + tau * ( -a*v - nablaf(x(:,end))  )/2;
            x(:,end+1) = x(:,end) + tau*v;
            v = v + tau * ( -a*v - nablaf(x(:,end))  )/2;
    end
    
end


q = 100; % #display
disp_list = round(linspace(1,niter,q));
for it=1:q
    %
    j = disp_list(it);
    clf; hold on;
    imagesc(t,t,rescale(f(Z))');
    contour(t,t,rescale(f(Z))',linspace(0,1,r), 'k');
    for k=1:K
        plot(x(k,1:j), 'color', col(k,:), 'LineWidth', 2);
        plot(real(x0(k)), imag(x0(k)), '.', 'color', col(k,:), 'MarkerSize', 15);
    end
    colormap(CM);
    caxis([0 1]);
    axis image; axis off;
    drawnow;
    mysaveas(it);
end



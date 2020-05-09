%% 
% Variational version of SOM.

n = 3000; % data points

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

addpath('../displ-interp-2d/toolbox-lsap/');
addpath('../displ-interp-2d/mexEMD/'); 

% generate data
c = [.7+.5i, .2+.3i, .15+.7i];
s = [.2, .1, .08]/2;
%
c = [.5+.5i];
s = [.1];
%
c = c(:); s = s(:);
%
I = ceil( rand(n,1)*length(c) );
x = randn(n,1)+1i*randn(n,1);
%x = (2*rand(n,1)-1)+1i*(2*rand(n,1)-1);
x = x.*s(I(:)) + c(I(:));

% grid init
m = 20; 
t = linspace(0,1,m);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;

I = randperm(n); 
%Z = reshape(x(I(1:m*m)), [m m]);

tau = 1; % increase for more speed evolution
lambda = .4; % decrease for more smooth

q = 50;
for it=1:q
    if 1 % mod(it,10)==1
    % plot grid
    clf; hold on;
    plot(x, '.');
    % plot(c, 'k.','MarkerSize', 25);
    %
    u = Z(1:end-1,:); v = Z(2:end,:);
    plot( real([u(:) v(:)])', imag([u(:) v(:)])', 'r-' );
    u = Z(:,1:end-1); v = Z(:,2:end);
    plot( real([u(:) v(:)])', imag([u(:) v(:)])', 'r-' );
    axis equal; axis([0 1 0 1]); axis off;
    drawnow;
    mysaveas(it);
    end
    
    % nn match
    D = abs(x(:)-transpose(Z(:)));  
	[cost,gamma] = mexEMD(ones(n,1)/n,ones(m*m,1)/(m*m),D.^2);
	g = reshape(gamma' * x, [m m]) * m*m;
    %
    t = .3;
    Z = .93 * Z + .07 * ( t*g + (1-t)*mysmoothing(g, 1) );
    
    % replace by an OT map
    
    
    if 0
    % evolve
    Znew = (1-tau)*Z + tau*g;
    Z = lambda*Znew + (1-lambda)*H;  
    end
end

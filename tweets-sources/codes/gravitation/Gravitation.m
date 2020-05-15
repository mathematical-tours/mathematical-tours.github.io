%%
% Display the core of a polygon

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

R = 1;

% initial speed and position




% force
F = @(x)-x./abs(x).^2;
v0 = (1+7*1i);
v0 = 3*v0/abs(v0)*[.3 .4 .45 .5 .53 .55]';
K = length(v0);
x0 = 1.2*R*ones(K,1);
T = 50;

% force -- true gravitation
F = @(x)-x./abs(x).^3;
v0 = (1+7*1i);
v0 = 2*v0/abs(v0)*[.3 .4 .45 .5 .54 .57]';
K = length(v0);
x0 = 1.2*R*ones(K,1);
T = 35;

% force
F = @(x)-x./abs(x).^4;
v0 = (1+7*1i);
v0 = 2*v0/abs(v0)*[.3 .4 .41 .415 .416 .42]';
K = length(v0);
x0 = 1.2*R*ones(K,1);
T = 100;




% force 
F = @(x)-x./abs(x);
v0 = (1+4*1i);
v0 = 5*v0/abs(v0)*[.1 .13 .2 .3 .4 .5]';
K = length(v0);
x0 = 1.2*R*ones(K,1);
T = 40;

% force harmonic
F = @(x)-x;
v0 = (1+4*1i);
v0 = 10*v0/abs(v0)*[.1 .13 .2 .3 .4 .5]';
K = length(v0);
x0 = 1.2*R*ones(K,1);
T = 10;

% 
tau = .01;
niter = round(T/tau);

x = x0; v = v0;
for j=1:niter
    x(:,end+1) = x(:,end) + tau/2 * v(:,end);
    v(:,end+1) = v(:,end) + tau * F(x(:,end));
    x(:,end) = x(:,end) + tau/2 * v(:,end);
end
for k=1:K
    I = find(abs(x(k,:))<R);
    if not(isempty(I))
        x(k,I(1):end) = NaN; % remove points collapsing
    end
end


theta = linspace(0,2*pi,100);
col = distinguishable_colors(K+1);
col(4,:) = [];
%
q = 200; % #display
disp_list = round(linspace(1,niter,q));
for it=1:q
    %
    j = disp_list(it);
    clf; hold on;
    fill( R*cos(theta),R*sin(theta), [1 1 1]*.2 );
    % contour(xlist,ylist,rescale(V(X,Y))',linspace(0,1,r), 'k');
    for k=1:K
        plot(real(x(k,1:j)), imag(x(k,1:j)), 'color', col(k,:), 'LineWidth', 2);
    end
    plot(real(x0), imag(x0), '.', 'color', 'k', 'MarkerSize', 20);
    axis equal; 
    axis([-1 1 -1 1]*5);  box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end




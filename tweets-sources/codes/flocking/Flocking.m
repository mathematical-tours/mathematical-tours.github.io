%%
% Cucker-Smale flocking model

addpath('../toolbox/');
rep = MkResRep();

% #particules
n = 100; 

% initial position and velocities
crandn = @(k)randn(k,1)+1i*randn(k,1);
x0 = crandn(n);
v0 = crandn(n) * .06+.02;

% x0 = [crandn(n/2)*.2 - 1.5;crandn(n/2)*.2 + 1.5];
% v0 = crandn(n) * .1;

r0 = .7;
psi = @(r)1./(1+(r/r0).^2).^.6 ;

% distance matrix
D = @(x)abs( repmat(x,[1 length(x)]) - repmat(x',[length(x) 1]) );
% laplacian
L = @(x)1/length(x)*( psi(D(x))-diag(sum(psi(D(x)))) );

% evolution in time
q = 100; % #frame
tau = 1.3; % step size
x = x0; v = v0; 
K = .2;
for i=1:q
    %
    c = (i-1)/(q-1);
    %
    [x,v] = deal(x+tau*(v), v + tau*K*L(x)*v); 
    % display
    x1 = ( x-mean(x) ) / std(x);
    x1 = x;
    s = v/median(abs(v))*.6;
    clf; hold on;
    plot( x1, 'k.', 'MarkerSize', 25, 'color', [c 0 1-c]);
    plot( transpose([x1 x1+s]), '-', 'LineWidth', 1, 'color', [c 0 1-c]);
    axis equal; axis off;
    axis([-1 1 -1 1]*6);
    drawnow;     
    % kill away particles
    % x = x(abs(x1)<=2.5); 
    % v = v(abs(x1)<=2.5);
    saveas(gcf, [rep 'flocking-' znum2str(i,3) '.png']);
end

% AutoCrop(rep, ['flocking-']); 
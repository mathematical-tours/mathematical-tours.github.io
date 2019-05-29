%%
% Heat equation on polynomials.

if not(exist('test'))
    test=1;
end
addpath('../toolbox/');
rep = MkResRep(num2str(test));


% second order differential on polynomials.
D2 = @(N)diag((1:N-2).*(2:N-1),+2); 

% initial data
P0 = [1 0 -1]'; % X^2-1;
P0 = [1 0 0 0 0 -1]'; % X^5-1;
P0 = [1 0 -1]'; % X^2-1;
P0 = [1 0 0 -1]'; % X^3-1;
P0 = [1 0 0 0 -1]'; % X^3-1;

%%% from roots

% random
v = 2*rand(6,1)-1;
tmax = .5;
% click and play
clf; hold on;
v = [];
while true
    axis equal; axis([-1 1 -1 1]); 
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    if button==3
        break;
    end
    v(end+1) = a+1i*b;
end
tmax = .5;

x= sym('x');
P0 = double( coeffs(expand(prod(x-v)))' );
N = length(P0);

myplot = @(r,m,col)plot(real(r), imag(r), '.', 'MarkerSize', m, 'Color', col);


% use simple explicit euler scheme.
niter = 500; 
tau = tmax/niter;
P = P0; R = [];
for i=1:niter
    t = (i-1)/(niter-1);
    P = P + tau*D2(N)*P;
    R(:,end+1) = roots(P(end:-1:1));
end

% grid for polyevaluation
Bx = 1; By = 2; m = 200;
xgrid = linspace(-Bx,Bx,m);
ygrid = linspace(-By,By,m);
[Y,X] = meshgrid(ygrid,xgrid); XY = X+1i*Y;


% use simple explicit euler scheme.
ms = 15; 
q = 50; 
ndisp = round(linspace(1,niter,q));
for it=1:q
    i = ndisp(it);
    clf; hold on;
    % polynomial background
    F = mypolyval(XY,R(:,i));
    U = log(.001+abs(F));
    U = max(U,-3);
    % display quantized colormap
    r = 15; % #levellines
    clf; hold on;
    imagesc(xgrid,ygrid,U');
    contour(xgrid,ygrid,U',linspace(min(U(:)),max(U(:)),r), 'k');
    %
    t = (it-1)/(q-1);
    mi = linspace(0,1,r-1)';
    CM = mi*[1 1 1] + (1-mi)*[t 0 1-t];
    %
    colormap(CM);    % parula(r-1) 
    caxis([min(U(:)) max(U(:))]);
    axis image; axis off;
    % plot(R(:,i), 'r.', 'MarkerSize', 25);
    %
    for j=1:i
        t = (j-1)/(niter-1);
        myplot( R(:,j), ms, [t 0 1-t]);
    end
    axis([-Bx Bx -By By]); axis equal; box on;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end

test = test+1;

% AutoCrop(rep, 'anim')
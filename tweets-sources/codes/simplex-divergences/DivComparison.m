%%
% Compares various divergences on 3D simplex.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep  'a-' name '-' znum2str(it,3) '.png']);

n = 256*2;
x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);
% (1,0)->(1,0), (0,1)->(cos(pi/3), sin(pi/3)
U = [1, cos(pi/3); 0, sin(pi/3)];
V = inv(U);
S1 = V(1,1)*X+V(1,2)*Y;
S2 = V(2,1)*X+V(2,2)*Y;
S3 = 1-S1-S2;
S = [S1(:) S2(:) S3(:)];
tau = 1e-3;
Sigma = (S1>=tau & S2>=tau & S3>=tau);
KL = @(x,y)sum(x.*log((x+1e-15)./(y+1e-15)),2);
Hell = @(x,y)sum((sqrt(x)-sqrt(y)).^2,2);
TV = @(x,y)sum(abs(x-y),2); vmax = 2;
PhiDiv = @(phi,x,y)abs( sum( phi(x./(y+1e-20)).*y ,2) );

if not(exist('name'))
name = 'tv';
name = 'hell';
name = 'kl';
name = 'klrev';
name = 'chi2';
name = 'js';
end

phi = [];
h = @(s)s; % post-scaling
switch name
    case 'tv'
        f = @(x,y)TV(x,y); vmax=2;
        phi = @(r)abs(1-r);
    case 'hell'
        f = @(x,y)sqrt(Hell(x,y)); vmax=sqrt(2);
        h = @(s)sqrt(s);
    case 'kl'
        f = @(x,y)sqrt(KL(x,y)); vmax=2.5;
        phi = @(r)r.*log(r+1e-15);
        h = @(s)sqrt(s);
    case 'klrev'
        f = @(x,y)sqrt(KL(y,x)); vmax=2.5;
        h = @(s)sqrt(s);
    case 'chi2'
        phi = @(r)(1-r).^2;
        h = @(s)sqrt(s);
        vmax = 3;
    case 'js'
        f = @(x,y)sqrt(KL(x,(x+y)/2) + KL(y,(x+y)/2) ); vmax=f([1 0 0],[0 1 0]);
        h = @(s)sqrt(s);
end
if not(isempty(phi))
    f = @(x,y)h(PhiDiv(phi,x,y)); 
end

% path
q = 100;
P = [1 1 1; .5 .5 0; 1 0 0; 1 1 0; 0 1 1; 0 1 0];
P = [1 1 1; 1 1 0; 1 0 1; 0 1 1; 0 1 0; 1 1 1];
P = P./sum(P,2);
Pint = interp1(linspace(0,1,size(P,1)),P,linspace(0,1,q), 'linear');
%
Pint = interp1(linspace(0,1,size(P,1)),P,linspace(0,1,q), 'spline');
Pint = max(Pint,1e-4); Pint = Pint./sum(Pint,2);

r = 15; % #levellines
CM = parula(r);
for it=1:q
    s=Pint(it,:);
    A = reshape(abs(f( s, S )), [n n]);
    A(Sigma==1) = min(A(Sigma==1),vmax)/vmax; % for KL
    % turn into color image
    I = 1+floor(clamp(A)*(r-1));
    U = reshape(CM(I,:), [n n 3]);
    for i=1:3
        Us = U(:,:,i); Us(Sigma==0) = 1; U(:,:,i) = Us;
    end
    % display quantized colormap
    A(Sigma==0) = NaN;
    clf; hold on;
    imagesc(x,x,permute(U, [2 1 3]));
    contour(x,x,A',linspace(0,1,r), 'k');
    plot( s(1) + s(2)*exp(1i*pi/3), 'r.', 'MarkerSize', 25 );
    axis image; axis off;
    drawnow;
    mysaveas(name, it);
end

return;

% TV<=sqrt(2*KL)

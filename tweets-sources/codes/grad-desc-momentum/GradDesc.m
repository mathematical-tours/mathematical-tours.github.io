%%
% Gradient descent on a quadratic function.

rep = 'results/';
[~,~] = mkdir(rep);

% x^2+a*y^2
a = 2;

gmode = 'nesterov1';
gmode = 'nesterov';
gmode = 'heavyball';

% anisotropy
a = 12;
% initial point
z0 = .9 + 1i * .3;

% display
N = 256;
tx = linspace(-.4,.92,N);
ty = linspace(-.25,.32,N);
[Y,X] = meshgrid(ty,tx);
F = ( X.^2+a*Y.^2 )/2;
                

clf; hold on;
imagesc(tx,ty,-F');
contour(tx,ty,F', 15, 'k');

% descent step size
tau = .01/10;
niter = round(2/tau);
% momentum
q = 9;
switch gmode
    case  'heavyball'
        mulist = linspace(0,.98,q);
    case {'nesterov' 'nesterov1'}
        mulist = linspace(0,.99,q);
end
% mulist = 0*[1 1];



% gradient of F
gF = @(z) real(z) + 1i * a*imag(z);

lgd = {''};
Zsvg = [];
for imu=1:length(mulist)
    mu = mulist(imu);
    
    c = (imu-1)/(length(mulist)-1);
    
    
    % initial point
    z = z0;
    Z = [];
    v = 0;
    s = z; zold = z;
    %v = -tau*gF(z);
    for i=1:niter
        Z(end+1) = z;
        % gradient
        switch gmode
            case  'heavyball'
                g = gF(z);
                v = mu*v - tau*g;
                z = z + v;
            case 'nesterov'
                g = gF(z + mu*v);
                v = mu*v - tau*g;
                z = z + v;
        end
    end
    Zsvg(:,end+1)=Z(:);
    
    plot(Z, '-', 'LineWidth', 2, 'MarkerSize', 23, 'Color', c*[1 0 0] + (1-c)*[0 0 1]);
    colormap parula(256);
    % axis image; 
    axis off; axis equal;
    axis([min(tx) max(tx) min(ty) max(ty)]);
    drawnow;
    lgd{end+1} = ['\mu=' num2str(mu)];
end
plot(0,0, 'k.', 'MarkerSize', 23);

saveas(gcf, [rep 'gd-' gmode  '.png'], 'png');
% legend(lgd);

if 0
clf;
plot(real(Zsvg));
end
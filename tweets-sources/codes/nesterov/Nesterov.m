%%
% Gradient descent on a quadratic function.

rep = '../results/nesterov/';
[~,~] = mkdir(rep);

% x^2+a*y^2
% anisotropy
a = 12;
% initial point
z0 = .9 + 1i * .3;


% anisotropy
a = 6;
% initial point
z0 = .3 + 1i * .3;

% display
N = 256;
tx = linspace(-.1,real(z0)*1.1,N);
ty = linspace(-.1,imag(z0)*1.1,N);

[Y,X] = meshgrid(ty,tx);

F = ( X.^2+a*Y.^2 )/2;
                

% gradient descent step size
tau_g = .01/10;
niter = round(5/tau_g);
% Nesterov descent step size
tau_n = .01/20;


% gradient of F
gF = @(z) real(z) + 1i * a*imag(z);



% initial point
z_g = z0;
z_n = z0; y_n = z0;
Z_g = [];
Z_n = []; 
%v = -tau*gF(z);
for k=1:niter
    Z_g(end+1) = z_g;
    Z_n(end+1) = z_n;
    % gradient
    z_g = z_g - tau_g * gF(z_g);
    % nesterov
    z_old = z_n;
    z_n = y_n - tau_n * gF(y_n);
    y_n = z_n + k/(k+3) * (z_n-z_old);
end

clf; hold on;
imagesc(tx,ty,-F');
contour(tx,ty,F', 15, 'k');
%
plot(Z_g, 'b-', 'LineWidth', 2);
plot(Z_n, 'r-', 'LineWidth', 2);
colormap parula(256);
% axis image;
axis off; axis equal;
axis([min(tx) max(tx) min(ty) max(ty)]);
drawnow;

plot(0,0, 'k.', 'MarkerSize', 23);

saveas(gcf, [rep 'nesterov.png'], 'png');
% legend(lgd);

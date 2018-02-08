%%
% Gradient descent on a quadratic function.

rep = 'results/';
[~,~] = mkdir(rep);

% x^2+a*y^2
a = 2;

gmode = 'large';
gmode = 'search';
gmode = 'low';
niter = 10;


for a=4 % [1 2 4 8 16]

N = 256;
tx = linspace(-.3,1,N);
ty = linspace(-.6,.6,N);
[Y,X] = meshgrid(ty,tx);

F = ( X.^2+a*Y.^2 )/2;


% initial point
x = .9; y = .3;
U = []; V = [];
for i=1:niter
    U(end+1) = x; V(end+1) = y;
    % min_t (x-t*gx)^2+a*(y-t*gy)^2
    %    (x-t*x)^2+a*(y-t*a*y)^2
    %    x^2*(1-t)^2+a*y^2*(1-a*t)^2
    %    x^2*(1-t)+y^2*a^2*(1-a*t)=0
    %   t*(x^2+a^3*y^2) = x^2+a^2*y^2
    %   t = (x^2+a^2*y^2)/(x^2+a^3*y^2)
    switch gmode
        case 'search'
            r = (x^2+a^2*y^2)/(x^2+a^3*y^2);
        case 'low'
            r = .1;
        case 'large'
            r = .52;
    end
    x = x - r*x; 
    y = y - r*a*y;
end

clf; hold on;
imagesc(tx,ty,F');
contour(tx,ty,F', 15, 'k');
plot(U,V, 'r.-', 'LineWidth', 2, 'MarkerSize', 23);
colormap parula(256);
axis image; axis off;
axis([min(tx) max(tx) min(ty) max(ty)]);
drawnow;


saveas(gcf, [rep 'gd-' num2str(a) '-' gmode '.png'], 'png');

end

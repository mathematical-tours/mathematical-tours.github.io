%%
% Display of interpolation between two Gaussians.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


R = @(t)[cos(t),sin(t);-sin(t),cos(t)];
C = @(t,a)R(t)*diag([1,a])*R(-t);

A = C(pi/3,.02); B = C(2*pi/3,.05);

A = C(pi/3,.01); B = C(pi/3+pi/2,.01);


n = 201;
t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
d = @(A)reshape( sum( [X(:)';Y(:)'] .* (inv(A)*[X(:)';Y(:)']) ), [n n]);
s = .4;
G = @(A)exp(-d(A)/(2*s^2));

r = 12; % #colors
m = linspace(0,1,r-1)';


q = 50; 
for it=1:q    
    s = (it-1)/(q-1);
    % interpolating
    U = inv(sqrtm(A)) * sqrtm( sqrtm(A)*B*sqrtm(A) ) * inv(sqrtm(A));
    H = ((1-s)*eye(2)+s*U) * A * ((1-s)*eye(2)+s*U);
    %
    clf; hold on;
    imagesc(t,t,G(H));
    contour(t,t,G(H),linspace(0,1,r), 'k');
    colormap( m*[s 0 1-s] + (1-m)*[1 1 1] );
    caxis([0 1]);
    axis equal; box on; set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('anim',it);
end

% AutoCrop(rep, 'anim-'); 


addpath('../toolbox/');
addpath('./toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);


J = 7;
options.relaxation = 2;
[V,F] = compute_semiregular_sphere(J,options);

% in spherical coords
[th,phi,R] = cart2sph(V(1,:),V(2,:),V(3,:));
n = size(V,2);

%% Display all harmonics

L = 10; % frequency level
z = .8;
% compute spherical harmonics
y = [];
for M=-L:L
    y(end+1,:) = sphrharm(phi,th,L,M);
end

y = real(y)';

q = 100;
p = 5; % # steps
w = randn(p,2*L+1); w = w./sqrt(sum(w.^2,2));
W = interp1(linspace(0,1,p),w,linspace(0,1,q), 'spline');

for it=1:q
    t = (it-1)/(q-1);
    z = sum( W(it,:) .* y,2  );
    clf;
    sphere_disp(V,F,real(z), real(z), 1, 1);
    view(3); camlight;
    axis equal; axis([-1 1 -1 1 -1 1]*2);
    colormap parula(256);
    set(gca, 'CameraViewAngle', 10)
    drawnow;
    mysaveas(it);
end
% saveas(gcf, [rep 'spharm-' num2str(L) '-' num2str(M) '.png' ]);

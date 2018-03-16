rep = '../results/spherical-harmonics/';
[~,~] = mkdir(rep);
addpath('./toolbox/');
addpath('../toolbox/');

J = 7;
options.relaxation = 2;
[V,F] = compute_semiregular_sphere(J,options);

% in spherical coords
[th,phi,R] = cart2sph(V(1,:),V(2,:),V(3,:));
n = size(V,2);

%% Display all harmonics

Lmax = 4;
z = .8;
% compute spherical harmonics
for L=0:Lmax
    for M=-L:L
        y = sphrharm(phi,th,L,M);   
        clf;
        sphere_disp(V,F,real(y), real(y), 1, 1);    
        view(3); camlight; 
        set(gca, 'CameraViewAngle', 10)
        drawnow;
        saveas(gcf, [rep 'spharm-' num2str(L) '-' num2str(M) '.png' ]);
    end
end

return;

Lmax = 4;
% compute spherical harmonics
Y = [];
for L=0:Lmax
    for M=-L:L
        Y(:,end+1) = sphrharm(phi,th,L,M);   
    end
end

view(3); camlight;


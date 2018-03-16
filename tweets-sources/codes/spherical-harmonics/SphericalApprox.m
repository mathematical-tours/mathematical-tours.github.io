rep = '../results/spherical-harmonics/approx/';
[~,~] = mkdir(rep);
addpath('./toolbox/');
addpath('../toolbox/');

J = 7;
options.relaxation = 2;
[V,F] = compute_semiregular_sphere(J,options);

% in spherical coords
[th,phi,R] = cart2sph(V(1,:),V(2,:),V(3,:));
n = size(V,2);


%% Compute all harmonics

Lmax = 70;
Y = [];
for L=0:Lmax
    progressbar(L+1,Lmax+1);
    for M=-L:L
        Y(:,end+1) = sphrharm(phi,th,L,M);   
    end
end
imagesc(Y'*Y);
% orthogonalize

[Ua,S,Va] = svd(Y, 'econ');
Y1 = Ua*Va';


%% Non-linear approximation
name = 'map-bw';
m = 256;
f0 = load_image(name, m);
f0 = rescale(mean(f0,3))';
f0 = f0(:,end:-1:1);
% resample
t = linspace(-pi,pi,m); p = linspace(-pi/2,pi/2,m);
[P,T] = meshgrid(p,t);
f = interp2(P,T,f0,phi,th);
f = f(:);

% initial
rho = .2;
clf; 
sphere_disp(V,F,f,f,1,rho);
view(73,-26); camlight; colormap jet(256);
saveas(gcf, [rep 'original.png' ]);


% approximation for several L

% number atoms per L
NbL = @(L)2 * L*(L+1)/2 + (L+1);

for L1 = 10:10:Lmax
    p = NbL(L1);
    fP = Y1(:,1:p) * (Y1(:,1:p)'*f);
    clf;
    sphere_disp(V,F,fP,fP,1,rho);
    view(73,-26); camlight; colormap jet(256);
    drawnow;
    saveas(gcf, [rep 'approx-' num2str(L1) '.png' ]);
end

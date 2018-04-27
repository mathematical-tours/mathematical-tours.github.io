addpath('./toolbox/');
rep = MkResRep('approx');

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
name = 'disks';
m = 256;
switch name
    case 'map-bw'        
        f0 = load_image(name, m);
        f0 = rescale(mean(f0,3))';
        f0 = f0(:,end:-1:1);
        % resample
        t = linspace(-pi,pi,m); p = linspace(-pi/2,pi/2,m);
        [P,T] = meshgrid(p,t);
        f = interp2(P,T,f0,phi,th);
        f = f(:);
    case 'disks'
        % use V positions
        r = .4;
        f = (sum(V(1:2,:).^2)<=r^2)  + ...
            (sum(V(2:3,:).^2)<=r^2) + ...
            (sum(V([1 3],:).^2)<=r^2);
        f = f(:);
        
end

% initial
rho = .2;
clf; 
sphere_disp(V,F,f,f,1,rho);
view(73,-26); camlight; colormap jet(256);
saveas(gcf, [rep name '-original.png' ]);


% approximation for several L

% number atoms per L
NbL = @(L)2 * L*(L+1)/2 + (L+1);

for L1 = 1:min(50,Lmax) % 10:10:Lmax
    p = NbL(L1);
    fP = Y1(:,1:p) * (Y1(:,1:p)'*f);
    clf;
    sphere_disp(V,F,fP,fP,1,rho);
    view(73,-26); camlight; colormap jet(256);
    drawnow;
    saveas(gcf, [rep name '-' 'approx-' znum2str(L1,2) '.png' ]);
end


% AutoCrop(rep, [name '-' 'approx-'])

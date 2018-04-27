%%
% Display wave/heat on meshes.

edp_name = 'wave';
edp_name = 'heat';

addpath('../toolbox/');
rep = MkResRep(edp_name);


name = 'elephant';
[X0,F] = read_off([name '.off']);
n = size(X0,2);
clear options; options.name = name;


clf;
plot_mesh(X0,F,options);
camlight;

[L,W,D,iD] = MeshLaplacian(X0,F);
tW = iD * W;
tL = iD * L;

% Generate Data on the surface
n = size(X0,2);
% width of the Gaussian
sigma = .01;
% numer of points
npoints = 16;
% initial condition
f0 = zeros(n,1);
% previous distances
D = zeros(n,1) + 1e10;
for i=1:npoints
    % select a point farthest away
    [tmp,p] = max(D);
    % compute the distance beetween the point and the other vertices
    d = sum( (X0 - X0(:,p)*ones(1,n)).^2 );
    f0 = f0 + (-1)^mod(i,2) *exp( -d(:)/(2*sigma.^2) );
    D = min(D,d(:));
end

% display
clf;
options.face_vertex_color = f0(:);
plot_mesh(X0,F, options);
camlight; colormap jet(256);

switch edp_name
    case 'heat'
        Tmax = 800;
        tau = .9;
        rvmax = .4; % memory for smoothing the clamping
        rho_clamp = .5;
    case 'wave'
        Tmax = 100;
        tau = .01/20; 
        rvmax = .4; % memory for smoothing the clamping
        rho_clamp = .5;
end
niter = round(Tmax/tau);
ndisp = round(niter/140);

f1 = f0;
f = f0; % init
vmax = max(abs(f)); 
clf;
for i=1:niter    
    switch edp_name
        case 'heat'
            f = f - tau*tL*f;
        case 'wave'
            [f,f1] = deal(2*f - f1 - tau^2 * tL*f, f);
            
    end
        vmax =  (1-rvmax) * vmax + rvmax * max(abs(f));
    if  ndisp==1 || mod(i,ndisp)==1
        g = f; 
        g(g==max(g(:))) = vmax;
        g(g==min(g(:))) =-vmax;
        g = clamp(g,-rho_clamp*vmax, rho_clamp*vmax);
        clf;
        options.face_vertex_color = g(:);
        plot_mesh(X0,F, options);
        camlight; colormap jet(256); zoom(.95);
        drawnow;
        saveas(gcf, [rep edp_name '-' znum2str(1+(i-1)/ndisp,3) '.png']);
    end
end


% AutoCrop(rep, [edp_name '-']); 

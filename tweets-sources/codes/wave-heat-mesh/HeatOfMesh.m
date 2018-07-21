%%
% Heat equation diffusion of the mesh itself.

addpath('../toolbox/');
rep = MkResRep('heat-on-mesh');


name = 'elephant';
name = 'bunny'; 
[X0,F] = read_off([name '.off']);
n = size(X0,2);
clear options; options.name = name;


clf;
plot_mesh(X0,F,options);
camlight;

[L,W,D,iD] = MeshLaplacian(X0,F);
tW = iD * W; tL = iD * L;

% Generate Data on the surface
n = size(X0,2);



% display
clf;
plot_mesh(X0,F, options);
camlight; colormap jet(256);


Tmax = 800;
tau = .9;

rvmax = .4; % memory for smoothing the clamping
rho_clamp = .5;


niter = round(Tmax/tau);


q = 70; % frames
displist = round(linspace(1,niter,q));
k = 1;

X = X0;
for i=1:niter    
    if i==displist(k)
        t = (k-1)/(q-1);
        clf;
        options.face_vertex_color = [t 0 1-t];
        plot_mesh(X,F, options);
        camlight; colormap jet(256); zoom(.95);
        drawnow;
        saveas(gcf, [rep  'heat-' znum2str(k,2) '.png']);
        k = k+1;
    end
    X = X - tau*(tL*X')';
        % update laplacian
       % [L,W,D,iD] = MeshLaplacian(X,F); L = real(L); W = real(W); D = real(D); iD = real(iD);
       % tW = iD * W; tL = iD * L;
    
end

% AutoCrop(rep, 'heat-');


%%
% Display wave/heat on meshes.

edp_name = 'wave';
edp_name = 'heat';

addpath('../toolbox/');
rep = MkResRep(edp_name);


name = 'elephant';
name = 'bunny';
[X0,F] = read_off([name '.off']);
n = size(X0,2);
clear options; options.name = name;


clf;
plot_mesh(X0,F,options);
camlight;

[L,W,D,iD] = MeshLaplacian(X0,F);
tW = iD * W;
tL = iD * L;
n = size(X0,2);

fname = 'lena';
fname = 'hibiscus';

% Generate Data on the surface
switch fname
    case 'bumps'
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
    otherwise
        p = 256;
        M0 = load_image(fname,p);
        M0 = rescale(sum(M0,3))';
        % equalize
        [fS,I] = sort(M0(:)); M0(I) = linspace(0,1,length(M0(:)));
        % zero pad
        M = zeros(round(1.4*p/2)*2) + 1/2;
        M(end/2-p/2+1:end/2+p/2,end/2-p/2+1:end/2+p/2) = M0;
        
        [Y,X] = meshgrid(1:p,1:p);
        M = double( (mod(X,18)<=8) | (mod(Y,18)<=8) );
        
        % sph coords
        v = X0 - repmat(mean(X0,2), [1 n]);
        theta = acos(v(1,:)./sqrt(sum(v.^2)))/pi;
        phi = (atan2(v(2,:),v(3,:))/pi+1)/2;
        % interp
        x = linspace(0,1,size(M,1));
        f0 = interp2(x,x,M',theta,phi)';
end

% display
clf;
options.face_vertex_color = f0(:);
plot_mesh(X0,F, options);
camlight; colormap gray(256);

switch edp_name
    case 'heat'
        Tmax = 100;
        tau = .2;
        rvmax = .4; % memory for smoothing the clamping
        rho_clamp = .5;
    case 'wave'
        Tmax = 100;
        tau = .01/20; 
        rvmax = .4; % memory for smoothing the clamping
        rho_clamp = .5;
end
niter = round(Tmax/tau);


q = 70; % frames
displist = round(linspace(1,niter,q));
k = 1;


f1 = f0;
f = f0; % init
vmax = max(abs(f)); 
clf;
for i=1:niter    
        vmax =  (1-rvmax) * vmax + rvmax * max(abs(f));
    if i==displist(k)
        g = f; 
        g(g==max(g(:))) = vmax;
        g(g==min(g(:))) =-vmax;
        g = clamp(g,-rho_clamp*vmax, rho_clamp*vmax);
        %
        g = rescale( clamp(rescale(f),.2,.8) );
        clf;
        options.face_vertex_color = g(:);
        plot_mesh(X0,F, options);
        camlight; colormap parula(256); zoom(.95);
        drawnow;        
        saveas(gcf, [rep  edp_name '-' znum2str(k,2) '.png']);
        k = k+1;        
    end
    switch edp_name
        case 'heat'
            f = f - tau*tL*f;
        case 'wave'
            [f,f1] = deal(2*f - f1 - tau^2 * tL*f, f);
            
    end
end


% AutoCrop(rep, [edp_name '-']); 

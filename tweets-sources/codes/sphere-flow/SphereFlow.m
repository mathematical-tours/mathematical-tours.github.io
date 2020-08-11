%%
% Display integration of a curve on a sphere.

addpath('../toolbox/');
addpath('../spherical-wavelets/toolbox_multires/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


% compute a multiresolution sphere.
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;
J = 5;
[VL,FL] = compute_semiregular_sphere(J,options);
V = VL{end}; F = FL{end};
Vh = VL{end}; Fh = FL{end};


% generate a random vector field in 3D
n = 50; 
x = [0:n/2,-n/2+1:-1]/n; x2=x.^2;
s = .16;
s = .04;
K = exp( -(x2+x2'+reshape(x2,[1 1 n]))/(2*s^2)  );
f = [];
for i=1:3
    f(:,:,:,i) = real(ifftn( fftn(randn(n,n,n)) .* fftn(K) ) );
end
vfield = @(V)Interp3D(f,V);

% variatinal setting (function minimization).
ProjSph = @(g,V)g - sum(g.*V) .* V;

vfield = @(V)ProjSph( V - .5*[-1 -1 -1]',V );
Func = @(V)sum( -( V - .5*[-1 -1 -1]' ).^2, 1);

v1 = [1 0 1]'; v1 = v1/norm(v1);
v2 = [-1 .5 4]'; v2 = v2/norm(v2);

f = @(V)sum((V-v1).^2); fD = @(V)V-v1;
g = @(V)sum((V-v2).^2); gD = @(V)V-v2;
normalize = @(g)g ./ sqrt(sum(g.^2)); 
vfield = @(V)-normalize( ProjSph( f(V).*gD(V) + g(V).*fD(V) , V ) );
Func = @(V)f(V).*g(V);

vfield = @(V)-normalize( ProjSph( gD(V) , V ) );
Func = @(V)g(V);


clf; hold on;
opt.face_vertex_color = ( rescale(Func(Vh))' ).^.2;
plot_mesh(Vh,Fh, opt);
%plot3(v1(1),v1(2),v1(3), 'r.', 'MarkerSize', 20);
% plot3(v2(1),v2(2),v2(3), 'r.', 'MarkerSize', 20);
colormap parula(256);
camlight;
saveas(gcf, [rep 'sphere-func.png']);

% interpolate field

if 0
g = vfield(V);
clf; 
DisplaySphereVF(Vh,Fh, V,g);
% saveas(gcf, [rep 'sphere-flow.png']);
end

% integrate the vf
q = 60; % #frames
tau = 1.5/q;
W = VL{4};
p = size(W,2);
col = distinguishable_colors(p);
for it=1:q
    % advance
    W(:,:,end+1) = W(:,:,end) + tau * vfield(W(:,:,end));
    % project on the sphere
    W(:,:,end) = W(:,:,end) ./ sqrt(sum(W(:,:,end).^2)); 
    % display
    clf; hold on;
    plot_mesh(Vh,Fh);
    % DisplaySphereVF(Vh,Fh, V,g);
    camlight;
    for i=1:p
        plot3(squeeze(W(1,i,:)),squeeze(W(2,i,:)),squeeze(W(3,i,:)), 'color', col(i,:), 'LineWidth', 2);
    end
    drawnow;
    mysaveas(it);
end

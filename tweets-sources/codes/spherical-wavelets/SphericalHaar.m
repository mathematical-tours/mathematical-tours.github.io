%%
% Test for the spherical haar transform.

addpath('../toolbox/');
addpath('./toolbox_multires/');
rep = MkResRep();


%%
% First compute a multiresolution sphere.

options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;
J = 8;
[vertex,face] = compute_semiregular_sphere(J,options);
n = size(face{end},2);

%%
% Display two examples of sphere.

for j=1:6
    clf;
    plot_mesh(vertex{j}, face{j});
    view(vv);
    shading faceted;
    camlight;
    drawnow;
    saveas(gcf, [rep 'sphere-' znum2str(j,2) '.png']);
end


%%
% Comput the center of each face.

x = [];
for i=1:3
    v = vertex{end}(i,:);
    x(i,:) = mean(v(face{end}));
end

name = 'hibiscus';
name = 'geometric';
n0 = 512;
M = rescale( sum(load_image(name, n0),3) );


%%
% Load a function on the sphere.
% Use the center of each face to sample the function.

f = rescale( load_spherical_function(M, x, options) );

%% 
% Display the function on the sphere.

vv = [170,70];
options.face_vertex_color = f;
clf;
plot_mesh(vertex{end}, face{end}, options);
view(vv);
colormap gray(256);
lighting none;


%% Compute the successive low pass approximations.

fj = f;
for j=1:6
    fj = reshape(fj, [length(fj)/4 4]);
    fj = mean(fj,2);
    clf;
    options.face_vertex_color = fj;
    plot_mesh(vertex{end-j}, face{end-j}, options);
    view(vv);
    colormap gray(256);
    shading faceted;
    camlight;
    drawnow;
    saveas(gcf, [rep name '-lowpass-' znum2str(j,2) '.png']);
end


%% Haar transform

[fw,nj] = PerformHaarSphere(f,+1,J);
f1 = PerformHaarSphere(fw,-1,J,nj);
norm(f-f1)/norm(f)

%% 
% Perform progressive Haar wavelet approximation.

q = 50; 
nlist = round(linspace(0.001,.2,q)*n);

u = .1*linspace(0,1,q).^2.5;
nlist = 1 + round(u*n);


fwS = sort(abs(fw(:)), 'descend');
for i=1:q
    n1 = nlist(i);
    T = fwS(n1+1);
    f1 = PerformHaarSphere( fw .* (abs(fw)>T) ,-1,J,nj );
    % display
    clf;
    options.face_vertex_color = clamp(f1);
    plot_mesh(vertex{end}, face{end}, options);
    view(vv);
    camlight;
    colormap gray(256);
    drawnow;
    saveas(gcf, [rep name '-approx-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, [name '-approx-']);

%% 
% Same with image haar.


Mw = perform_haar_transf(M, 1, +1);
MwS = sort(abs(Mw(:)), 'descend');
for i=1:q
    n1 = nlist(i);
    T = MwS(n1+1);
    M1 = perform_haar_transf(Mw .* (abs(Mw)>T), 1, -1);
    imwrite(clamp(M1), [rep name '-image-' znum2str(i,2) '.png']);
    % display
    clf;
    imageplot(M1);
    drawnow;
end
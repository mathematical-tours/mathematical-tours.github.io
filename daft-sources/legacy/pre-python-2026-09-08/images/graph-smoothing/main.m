% test of mesh parameterization
%
%   Copyright (c) 2004 Gabriel Peyré

global lang;
if ~strcmp(lang,'eng') % default is french
    lang = 'fr';
end

filename = '../mesh/mannequin.off';
filename = '../mesh/nefertiti.off';
% [vertex,face] = read_off(filename);
A = triangulation2adjacency(face);
xy = vertex(:,1:2);
nface = length(face);
nvert = max(max(face));
% ring = compute_1_ring(face);

clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test normal computations
normal = compute_normal(vertex,face);

% plot the normal
subplot(1,2,1);
plot_mesh(vertex,face,0.7,0,normal)
if strcmp(lang,'eng')==1
    title('Mesh with its normals');
else
    title('Le maillage et ses normales');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ajout de bruit
ns = 0.4;   % scale of the noise
n = (2*rand(nvert,1)-1)*ns;
vertex2 = vertex + [ normal(:,1).*n, normal(:,2).*n, normal(:,3).*n ];

% plot the mesh
subplot(1,2,2);
plot_mesh(vertex2,face)
if strcmp(lang,'eng')==1
    title('With noise added');
else
    title('Avec ajout de bruit');
end

if exist('save_image')
    save_image('graph-normal');
end

disp('Press any key to continue');
pause;

clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combinatorial Laplacian smoothing
plot_at = [1,2,5];
niter = max(plot_at);
lap = compute_laplacian(A);
epsilon = 1;
T = eye(nvert,nvert)-epsilon*lap;
k = 1;
% smooth the mesh
disp('--> Smoothing with epsilon=1.');
vertex3 = vertex2;
subplot(1,length(plot_at)+1,1);
plot_mesh(vertex3,face);
axis tight;
title('t=0');
h = waitbar(0,'Smoothing with epsilon=1.');
for i=1:niter
    waitbar(i/niter);
    vertex3 = T*vertex3;
    if i==plot_at(k)
        subplot(1,length(plot_at)+1,k+1);
        k = k+1;
        plot_mesh(vertex3,face);
        axis tight;
        title(sprintf('t=%d', i*epsilon));
    end
end
close(h);

if exist('save_image')
    save_image('graph-smoothing');
end

return

disp('Press any key to continue');
pause;



clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epsilon petit
plot_at = [1,2,5,10]*10;
niter = max(plot_at);
lap = compute_laplacian(A);
epsilon = 0.1;
T = eye(nvert,nvert)-epsilon*lap;
k = 1;
% smooth the mesh
disp('-->  with epsilon=0.1.');
vertex3 = vertex2;
h = waitbar(0,'Smoothing with epsilon=0.1.');
for i=1:niter
    waitbar(i/niter);
    vertex3 = T*vertex3;
    if i==plot_at(k)
        subplot(1,length(plot_at),k);
        k = k+1;
        plot_mesh(vertex3,face);
        title(sprintf('t=%.2f', i*epsilon));
        if k==2
            ylabel( sprintf('\epsilon=%.2f',epsilon) );
            set(get(gca,'YLabel'),'Visible','On');
        end
    end
end
close(h);


if exist('save_image')
    save_image('graph-smoothing-2');
end

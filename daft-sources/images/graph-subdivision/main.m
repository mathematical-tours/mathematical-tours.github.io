% Test of subdivision schemes
%
%   Copyright (c) 2003 Gabriel Peyré


global lang;
if ~strcmp(lang,'eng') % default is french
    lang = 'fr';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First show the subdivision process

mesh_types{1} = 'tetra';
mesh_types{2} = 'oct';
mesh_types{3} = 'ico';
nmesh = length(mesh_types);
nsub = 3;
clf;
set(gcf,'Renderer','zbuffer')
for i=1:nmesh
    mesh_type = mesh_types{i};
    mesh_name = mesh_type;
    [vertex,face] = gen_base_mesh(mesh_type);
    if i==1
        vertex = vertex*compute_rotation([0,1,0], -pi/2)';
        face = reverse_orientation(face);
    end
    if i==2
        vertex = vertex*compute_rotation([1,0,0], pi/8)';
        face = reverse_orientation(face);
    end
    for s=0:nsub
        subplot(nmesh,nsub+1,s+1+(nsub+1)*(i-1));
        plot_mesh(vertex,face);
        lighting flat;
        if i==1
            str = sprintf('Subdivision %d', s);
            title(str);
        end
        if s==0
            if strcmp(lang,'eng')==1
                ylabel( mesh_name );
            else
                ylabel( mesh_name );
            end
            set(get(gca,'YLabel'),'Visible','On');
        end
        if s~=nsub     
            [vertex,face] = subdivide_sphere(vertex,face);
        end
    end
end
    

if exist('save_image')
    save_image('graph-subdivision-sphere');
end

disp('Press a key to continue.');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a sphere
[vertex,face] = gen_base_mesh(mesh_type);
[vertex,face] = subdivide_sphere(vertex,face,3);
nvert = size(vertex,1);
nface = size(face,1);

% the adjacency matrix
A = triangulation2adjacency(face);

% build the laplacian from adjacency matrix
lap = compute_laplacian(A);

% performing SVD
[U,S,V] = svd(lap);
n = length(lap);

p = 3;
clf;
for i=1:p^2
    num = 5*(i-1)+1;
    c = U(:,n-num);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    subplot(p,p,i)
    
    patch('vertices',vertex,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    
    if strcmp(lang,'eng')==1
        str = ['Eigenvector n°', num2str(num+1)];
    else
        str = ['Vecteur propre n°', num2str(num+1)];
    end
    title(str);
end

if exist('save_image')
    save_image('graph-subdivision-vecteurs-propres');
end

return;
disp('Press a key to continue.');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Displaying reconstruction.');
% reconstruction 
p = 3;
nbr_max_keep = 20;
clf;
for i=1:p^2
    subplot(p,p,i)

    keep = round(i*nbr_max_keep/p^2); % nbr de pourcent gardé
    vertex2 = U'*vertex;
    % on détermine le seuil
    vnorm = vertex2(:,1).^2 + vertex2(:,2).^2 + vertex2(:,2).^2;
    vnorms = sort(vnorm); vnorms = reverse(vnorms);
    thresh = vnorms( round(keep/100*nvert) );
    % on enlève tous les coeff en dessous du seuil
    mask = vnorm>thresh;
    vertex2(:,1) = vertex2(:,1).*mask;
    vertex2(:,2) = vertex2(:,2).*mask;
    vertex2(:,3) = vertex2(:,3).*mask;
    vertex2 = U*vertex2;

    patch('vertices',vertex2,'faces',face,'facecolor',[cf cf cf],'edgecolor',[ce ce ce]);
    title([num2str(keep), '% des coefficients']);
    lighting phong;
    camlight infinite; 
    camproj('orthographic');
    axis square; axis off; 
end

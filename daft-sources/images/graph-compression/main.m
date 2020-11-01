% Teste la compression spectrale, tiré de l'article
%
%   Karni Z. and Gotsman C.
%   Spectral Compression of Mesh Geometry.
%   Computer Graphics (Proceedings of SIGGRAPH), pp. 279-286, 2000. 
%   <http://www.cs.technion.ac.il/~gotsman/AmendedPubl/SpectralCompression/SpectralCompression.pdf>
%
%   Copyright (c) 2004 Gabriel Peyré

global lang;
if ~strcmp(lang,'eng') % default is french
    lang = 'fr';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display some example of 3D meshes
base_path = '../mesh/';
clear filename;
filename{1} = [base_path, 'mushroom.off'];
filename{2} = [base_path, 'nefertiti.off'];
filename{3} = [base_path, 'venus.off'];
% filename{4} = [base_path, 'mannequin.off';
% filename{5} = [base_path, 'cube_oversampled.off';
% filename{6} = [base_path, 'pipes-loop2.off';
% filename{7} = [base_path, 'hhh.off';
nmesh = length(filename);

clf;
for i=1:nmesh
    subplot(1,nmesh,i);
    file = filename{i};
    [vertex,face] = read_off(file);
    plot_mesh(vertex,face);
    
    % retrieve name of the mesh
    k = strfind(file,'/');
    mesh_name = file((k(end)+1):end);
    k = strfind(mesh_name,'.');
    mesh_name = mesh_name(1:(k(1)-1));
    mesh_name = strrep(mesh_name,'_','\_');
    
    if strcmp(lang,'eng')==1
        str = ['Mesh ', mesh_name];
    else
        str = ['Maillage ', mesh_name];
    end
    title(str);
end

if exist('save_image')
    save_image('graph-examples-meshes');
end

disp('Press any key to continue.');
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same with normals
base_path = '../mesh/';
clear filename;
filename{1} = [base_path, 'mushroom.off'];
filename{2} = [base_path, 'nefertiti.off'];
filename{3} = [base_path, 'venus.off'];
% filename{4} = [base_path, 'mannequin.off';
% filename{5} = [base_path, 'cube_oversampled.off';
% filename{6} = [base_path, 'pipes-loop2.off';
% filename{7} = [base_path, 'hhh.off';
nmesh = length(filename);

clf;
for i=1:nmesh
    subplot(1,nmesh,i);
    file = filename{i};
    [vertex,face] = read_off(file);
    normal = compute_normal(vertex,face);
    plot_mesh(vertex,face,0.7,0,normal);
    
    % retrieve name of the mesh
    k = strfind(file,'/');
    mesh_name = file((k(end)+1):end);
    k = strfind(mesh_name,'.');
    mesh_name = mesh_name(1:(k(1)-1));
    mesh_name = strrep(mesh_name,'_','\_');
    
    if strcmp(lang,'eng')==1
        str = ['Mesh ', mesh_name];
    else
        str = ['Maillage ', mesh_name];
    end
    title(str);
end

if exist('save_image')
    save_image('graph-examples-meshes-normal');
end

disp('Press any key to continue.');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the 1st mesh
file = filename{1};
[vertex,face] = read_off(file);

nvert = size(vertex,1);
nface = size(face,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot eigenvectors
disp('Computing laplacian matrix.');
% the adjacency matrix
A = triangulation2adjacency(face);

% build the laplacian from adjacency matrix
lap = compute_laplacian(A);

disp('Performing SVD.');
tic;
[U,S,V] = svd(lap);
disp( sprintf('CPU time : %.2f.', toc) );

disp('Displaying eigenvectors.');
p = 3;
clf;
for i=1:p^2
    num = 14*(i-1)+1;
    c = U(:,nvert-num);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    subplot(p,p,i);
    patch('vertices',vertex,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    if strcmp(lang,'eng')==1
        str = ['Eigenvector n°', num2str(num+1)];
    else
        str = ['Vecteur propre n°', num2str(num+1)];
    end
    title(str);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
end


if exist('save_image')
    save_image('graph-compression-eigenvectors');
end

disp('Press any key to continue.');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot reconstruction
disp('Displaying reconstruction.');
% reconstruction 
p = 3;
nbr_max_keep = 8;
clf;
for i=1:p^2
    subplot(p,p,i);

    keep = 1+round(i*nbr_max_keep/p^2); % nbr de pourcent gardé
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

    plot_mesh(vertex2,face);
    if strcmp(lang,'eng')==1
        str = [num2str(keep), '% of the coefficients'];
    else
        str = [num2str(keep), '% des coefficients'];
    end
    title(str);
end

if exist('save_image')
    save_image('graph-compression-progression');
end
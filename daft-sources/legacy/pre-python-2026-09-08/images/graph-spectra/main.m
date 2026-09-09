% test for the spectra of some graph
%
%   Copyright (c) 2004 Gabriel Peyré

global lang;
if ~strcmp(lang,'eng') % default is french
    lang = 'fr';
end

if 0 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a triangle
nsub = 3;
[vertex,face] = gen_base_mesh('triangle');
clf;
for i=0:nsub
    subplot(1,nsub+1,i+1);
    plot_mesh(vertex,face);
    title( sprintf('Subdivision %d', nsub) );
    if i~=nsub
        [vertex,face] = subdivide(vertex,face);
    end
end

if exist('save_image')
    save_image('graph-subdivision-triangle');
end

disp('Press any key to continue.');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectra of a triangle
[vertex,face] = gen_base_mesh('triangle');
[vertex,face] = subdivide(vertex,face,5);
nface = size(face,1);
nvert = max(max(face));

% build the laplacian from adjacency matrix
disp('Computing Laplacian.');
A = triangulation2adjacency(face);
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
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    
    title(str);
end

if exist('save_image')
    save_image('graph-spectra-triangle');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the 30th eigenvectors for various shapes
shape_types{1} = 'square';
shape_types{2} = 'square1';
shape_types{3} = 'L1';
nsub = [3,4,2];
% nsub = [4,4,3];
num1 = 15;
num2 = 30;

nshape = length(shape_types);
clf;
for i=1:nshape
    shape_type = shape_types{i};
    
    [vertex,face] = gen_base_mesh(shape_type);
    [vertex,face] = subdivide(vertex,face,nsub(i));
    nface = size(face,1);
    nvert = max(max(face));
    
    % plot surface
    subplot(4,nshape,i);
    plot_mesh(vertex,face);
    axis tight;
    if i==1
        if strcmp(lang,'eng')==1
            ylabel( 'Manual drawing' );
        else
            ylabel( 'Dessin manuel' );
        end
        set(get(gca,'YLabel'),'Visible','On');
    end
    
    disp('Computing Laplacian.');
    A = triangulation2adjacency(face);
    lap = compute_laplacian(A);
    

    disp('Performing SVD.');
    tic;
    [U,S,V] = svd(lap);
    disp( sprintf('CPU time : %.2f.', toc) );
    
    % compute laplacian embeding
    vertex1 = U(:,(nvert-2):(nvert-1));
    vertex1 = [vertex1,zeros(nvert,1)];
    
    % plot laplacian embedding
    subplot(4,nshape,i+nshape);
    plot_mesh(vertex1,face);
    axis tight;
    if i==1
        if strcmp(lang,'eng')==1
            ylabel( 'Laplacian-based drawing' );
        else
            ylabel( 'Dessin à l''aide du Laplacien' );
        end
        set(get(gca,'YLabel'),'Visible','On');
    end
    
    % plot eigenvector n°1
    subplot(4,nshape,i+2*nshape);
    c = U(:,nvert-num1);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    patch('vertices',vertex1,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    axis tight;
    if i==1
        if strcmp(lang,'eng')==1
            ylabel( sprintf('Eigenvector n°%d',num1) );
        else
            ylabel( sprintf('Vecteur propre n°%d',num1) );
        end
        set(get(gca,'YLabel'),'Visible','On');
    end


    % plot eigenvector n°2
    subplot(4,nshape,i+3*nshape);
    c = U(:,nvert-num2);
    % rescale c
    c = (c-min(c))/(max(c)-min(c))*255;
    patch('vertices',vertex1,'faces',face,'FaceVertexCData',c,'edgecolor',[.2 .2 .6]);
    lighting phong; shading interp;
    colormap gray;
    axis square; axis off; 
    axis tight;
    if i==1
        if strcmp(lang,'eng')==1
            ylabel( sprintf('Eigenvector n°%d',num2) );
        else
            ylabel( sprintf('Vecteur propre n°%d',num2) );
        end
        set(get(gca,'YLabel'),'Visible','On');
    end
    
end

if exist('save_image')
    save_image('graph-spectra-diff-shape');
end
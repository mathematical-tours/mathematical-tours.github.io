% test of mesh parameterization
%
%   Copyright (c) 2004 Gabriel Peyré

global lang;
if ~strcmp(lang,'eng') % default is french
    lang = 'fr';
end

filename = '../mesh/mannequin.off';
filename = '../mesh/nefertiti.off';
[vertex,face] = read_off(filename);
A = triangulation2adjacency(face);
xy = vertex(:,1:2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute 1-ring
ring = compute_1_ring(face);


boundary_types{1} = 'circle';
boundary_types{2} = 'square';
boundary_types{3} = 'triangle';
nbound = length(boundary_types);
lap_types{1} = 'combinatorial';
lap_types{2} = 'conformal';
% lap_types{3} = 'authalic';
nlap = length(lap_types);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modèle d'origine
clf;
subplot(2,2,1);
cf = 0.8;
ce = 0;
patch('vertices',vertex,'faces',face,'facecolor',[cf cf cf],'edgecolor',[ce ce ce]);
lighting phong;
camlight infinite; 
camproj('perspective');
axis square; axis off;
title('Original model');
clf;

kk = 0;
for l = lap_types
k = 0;
kk = kk+1;
for b = boundary_types
    
    boundary_type = cell2mat(b);
    lap_type = cell2mat(l);
    
    k = k+1;
    str = sprintf('%s laplacian, %s boundary', lap_type, boundary_type);
    disp(['Computing parameterization : ', str, '.']);
    xy_spec = compute_parametrization(vertex,face,lap_type,boundary_type,ring);
    % essaie de redresse le graphe
    xy_spec = rectify_embedding(xy,xy_spec);

    subplot(1,nbound,k);
    gplot(A,xy_spec,'k.-');
    axis tight;
    axis square;
    axis off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up name for display
    if strcmp(lang,'eng')==1
        switch lap_type
        case 'combinatorial'
            lap_name = 'Combinatorial laplacian';
        case 'conformal'
            lap_name = 'Conformal laplacian';
        case 'authalic'
            lap_name = 'Authalic laplacian';
        end
        boundary_name = sprintf('%s boundary', boundary_type);
    else % french
        switch lap_type
        case 'combinatorial'
            lap_name = 'Laplacien combinatoire';
        case 'conformal'
            lap_name = 'Laplacien conforme';
        case 'authalic'
            lap_name = 'Laplacien authalique';
        end
        switch boundary_type
        case 'circle'
            boundary_name = 'Frontière circulaire';
        case 'square'
            boundary_name = 'Frontière carrée';
        case 'triangle'
            boundary_name = 'Frontière triangulaire';
        end
    end

    if kk==1
        title(boundary_name);
    end
    if mod(k-1,nbound)==0
        ylabel(lap_name);
        set(get(gca,'YLabel'),'Visible','On');
    end

end

if exist('save_image')
    save_image( sprintf('graph-parameterization-%d',kk) );
end

disp('Press any key to continue.');
pause;

end

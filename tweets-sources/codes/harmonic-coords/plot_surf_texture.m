function h = plot_surf_texture(M, T)

% plot_surf_texture - plot a surface with a texture on it.
%
%   h = plot_surf_texture(M, T);
%
%   M(:,:,i) for i=1,2,3 are the 3D coordonate of the surface.
%   T is a 2D image mapped on the surface.
%
%   Copyright (c) 2010 Gabriel Peyre

T = permute(T, [2 1 3]);
T = rescale(T);
if size(T,3)==1
    T = T*255;
end
h = surf(M(:,:,1), M(:,:,2), M(:,:,3) );
set(h,'CData',T,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct');
if size(T,3)==1
    colormap(gray(256));
else
   colormap(jet(256)); 
end
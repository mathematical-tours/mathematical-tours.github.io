function VolumetricRendering(f, color)


n = size(f,1);
x = linspace(-1,1,n);

F = isosurface( x,x,x,f, max(f(:))/2 );
p = patch(F);
isonormals( x,x,x,f,p );
set(p, 'FaceColor', color, 'EdgeColor', 'none');
alpha(1);
axis equal; axis([-1 1 -1 1]);
axis off; lighting phong;

end
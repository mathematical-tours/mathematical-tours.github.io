function VolumetricSlicing(f, slice_dim, color, slice_pos)

if nargin<2
    slice_dim = 1;
end
if nargin<4
    slice_pos=.5;
end

n = size(f,1);
x = linspace(-1,1,n);


R = 3; 
R = 8; % #slices
r = linspace(-1,1,R+2); r = r(2:end-1);
Col = jet(R);

for i=1:R
    t = (i-1)/(R-1);
    % mycolor = (1-r(i))*[1 0 1] + r(i)*[0 1 0];
    mycolor = (1-t)*color + t*[0 1 0];
    mycolor = Col(i,:);
    switch slice_dim
        case 1
            h = f; h(:,1:round(slice_pos*end),:) = NaN;
        case 2
            h = f; h(1:round(slice_pos*end),:,:) = NaN;
        case 3
            h = f; h(:,:,round(slice_pos*end)+1:end) = NaN;
    end
    F = isosurface( x,x,x,h, max(f(:))*r(i) );
    p = patch(F);
    isonormals( x,x,x,h,p );
    set(p, 'FaceColor', mycolor, 'EdgeColor', 'none');
end
alpha(1);
axis equal; axis([-1 1 -1 1]);
axis off; lighting phong;

end
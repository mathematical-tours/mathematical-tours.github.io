function PlotCube(C,I, col, lw, ms)

if nargin<3
    col = 'k';
end
if nargin<4
    lw = 2;
end
if nargin<5
    ms = 25;
end

Col = distinguishable_colors(size(I,2));

hold on;
for k=1:size(I,2)
    a = I(:,k);
    h = plot3( C(1,a), C(2,a), C(3,a), '-', 'LineWidth', lw, 'MarkerSize', ms, 'Color', Col(k,:) ); % 'Color', col
end
h = plot3( C(1,:), C(2,:), C(3,:), 'k.', 'MarkerSize', ms ); % 'Color', col
view(3);  axis off;


end
function RenderPolygon(Q,f,vmin,vmax,Z)

% display in color inside a polygon


p = size(f, 1);
t = linspace(0,1,p);
%
f = reshape(f,[p p]);
g = clamp(f,vmin,vmax);
g = rescale(g,0,1);

% render in colors the function
r = 20; % #levellines
c = parula(r);
I = floor((r-1)*g)+1;
h = c(I(:),:); 
h = reshape(h, [p p 3]);
% outside polygon
I = find(Q==0);
h([I,I+p^2,I+2*p^2])=1;
% display
hold on;
imagesc(t,t,permute(h,[2 1 3]));
plot(Z([1:end 1]), 'k.-', 'LineWidth', 2, 'MarkerSize', 25);
%
g1 = g; g1(Q==0) = NaN;
contour(t,t,g1',linspace(0,1,r), 'k');
%
axis equal; axis([0 1 0 1]); axis off;



end
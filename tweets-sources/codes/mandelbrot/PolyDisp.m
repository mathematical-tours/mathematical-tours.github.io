function PolyDisp(F,x,y)

n = size(F,1);
Fa = abs(F);

% levelset coords
r = 15; % #levellines
ls = linspace(0,max(Fa(:)),r);
%
% ls = sort(Fa(:)); ls = ls(round(linspace(1,n*n,r)));
% display quantized colormap
hold on;
imagesc(y,x,angle(F)');
% contour(x,y,abs(Fa)',ls, 'k');
colormap(hsv(256)); caxis([-pi pi]);
axis image; axis off;

end
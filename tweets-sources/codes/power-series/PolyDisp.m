function PolyDisp(F,t)

n = size(F,1);
Fa = abs(F);

% levelset coords
r = 15; % #levellines
%
ls = sort(Fa(:)); ls = ls(round(linspace(1,n*n,r)));
ls = linspace(0,max(Fa(:)),r);
% display quantized colormap
hold on;
imagesc(t,t,angle(F)');
contour(t,t,abs(Fa)',ls, 'k');
colormap(hsv(256)); caxis([-pi pi]);
axis image; axis off;

end
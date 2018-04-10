function DispOval(Z,z, t)

F = Z*0+1;
for i=1:length(z);
    F = F .* abs(Z-z(i));
end
F = F.^(1/length(z));

n = size(F,1);


% levelset coords
r = 15; % #levellines
% ls = sort(F(:)); ls = ls(round(linspace(1,n*n,r)));
% display quantized colormap
hold on;
imagesc(t,t,F');
contour(t,t,F',linspace(min(F(:)),max(F(:)),r), 'k');
colormap(parula(r-1));
caxis([min(F(:)),max(F(:))]);
% contour(t,t,F,ls, 'k');
% colormap(parula(256)); caxis([-pi pi]);
axis image; axis off;

end
%%
% Compare linear and non-linear approximation.

addpath('../toolbox/');


rep = ['../results/approximation/'];
[~,~] = mkdir(rep);

n = 512;

name = 'hibiscus';
f = load_image(name, n);
f = rescale(sum(f,3));

%%
% Linear approximation.

plist = [1:12, linspace(14,n,48)];
plist = 1:50;
for i=1:length(plist)
    f1 = LinApprox(f,plist(i));
    imwrite(f1, [rep 'lin-' num2str(i) '.png'], 'png');
    % flat
    clf; imageplot(f1); drawnow;
    % 3D
    clf; surf(f1);
    view(120,65); axis tight;
    shading interp; camlight; axis off; colormap gray(256); axis([1 n 1 n 0 1]);
    drawnow;
    % saveas(gcf, [rep '3d-lin-' num2str(i) '.png'], 'png');
end

%%
% Non linear approximation. 

fw = perform_haar_transf(f, 1, +1);
s = sort(abs(fw(:)), 'descend');
for i=1:length(plist)
    p = plist(i);
    fwT = fw .* (abs(fw)>=s(p*p));
    f1 = perform_haar_transf(fwT, 1, -1);
    % flat
    clf; imageplot(f1); drawnow;
    % imwrite(f1, [rep 'nlin-' num2str(i) '.png'], 'png');
    % 3D
    clf; surf(f1);
    view(120,65); axis tight;
    shading interp; camlight; axis off; colormap gray(256); axis([1 n 1 n 0 1]);
    drawnow;
    % saveas(gcf, [rep '3d-nlin-' num2str(i) '.png'], 'png');
end



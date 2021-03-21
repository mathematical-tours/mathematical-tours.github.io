

addpath('../toolbox/');
rep = MkResRep();

if not(exist('name'))
name = 'laplacian';
name = 'hilbert';
name = 'pascal';
name = 'gaussian';
end
switch name
    case 'laplacian'
        A = @(I,J)2*(I==J) - (I==J-1) - (I-1==J);
    case 'gaussian'
        s = 10; 
        A = @(I,J)exp( -(I-J).^2/(2*s));
    case 'hilbert'
        A = @(I,J)1./(I+J-1);
    case 'pascal'
        A = @(I,J)pascal(size(I,1));
end

q = 50;
for it=1:q    
    n = it+1;
    [J,I] = meshgrid(1:n,1:n);
    B = A(I,J);
    B = B/max(abs(B(:)));
    clf;
    imagesc(B);
    axis image; axis off;
    caxis([-1 1]);
    saveas(gcf, [rep 'matrix-' name '-' znum2str(it,3) '.png']);
    drawnow;
    
    a = real(eig(A(I,J)));
    clf;
    plot(1:n,a/max(a), 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
    axis([1 n -.03 1.03]);
    set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);    
    saveas(gcf, [rep 'eig-' name '-' znum2str(it,3) '.png']);
end


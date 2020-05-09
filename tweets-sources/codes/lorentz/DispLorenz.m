addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% RHO,SIGMA,BETA;
[x, y, z,T] = lorenz(28, 10, 8/3);
U = [x,y,z];
T = T/max(T);

U = U-mean(U);
U = U/max(abs(U(:)));

clf;
plot3(U(:,1),U(:,2),U(:,3));
n = size(U,1);
q = 150;

clf; hold on; 
plot3(U(:,1),U(:,2),U(:,3), 'color', [1 1 1]*0, 'LineWidth', .1);
for it=1:q
    s = it/q;
    I = find( (T>=(it-1)/q) & (T<=it/q) );
    I = [max(1,I(1)-1);I];
    plot3(U(I,1),U(I,2),U(I,3), 'color', [s 0 1-s], 'LineWidth', 3); 
    view(20,20); 
    axis tight; axis off;
    drawnow;
    mysaveas(it);
end
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

d = 3;
n = 20000;


v0 = 4/3*pi;

nruns = 500;
vx = zeros(n,nruns);
vy = zeros(n,nruns);
for j=1:nruns
    progressbar(j,nruns);
    p = sobolset(3);
    p = scramble(p,'MatousekAffineOwen');
    x = net(p,n);
    x = 2*x-1;
    y = 2*rand(n,d)-1;
    %
	vx(:,j) = 8*cumsum(sum(x.^2,2)<1) ./ (1:n)';
    vy(:,j) = 8*cumsum(sum(y.^2,2)<1) ./ (1:n)';
end

Ex = sqrt( mean( (vx-v0).^2, 2 ) );
Ey = sqrt( mean( (vy-v0).^2, 2 ) );


clf; 
semilogy((abs(vx(:,1)-v0)), 'LineWidth', 1, 'color', .5*[0 0 1] + .5);
hold on;
semilogy((abs(vy(:,1)-v0)), 'LineWidth', 1, 'color', .5*[1 0 0] + .5);
% mean
semilogy(Ex, 'LineWidth', 1, 'color', [0 0 1], 'LineWidth', 3);
semilogy(Ey, 'LineWidth', 1, 'color', [1 0 0], 'LineWidth', 3);
axis([0 n 1e-3 1]);
% axis([0 n -3 0]);
set(gca, 'XTick', [], 'YTick', [1e-3 1e-2 1e-1], 'FontSize', 13);
box on;
drawnow;
saveas(gcf, [rep 'mean-err.png']);
% mysaveas(it);


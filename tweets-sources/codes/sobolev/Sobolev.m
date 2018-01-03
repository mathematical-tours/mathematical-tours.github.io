%%
% Display random Sobolev functions

rep = '../results/sobolev/';
[~,~] = mkdir(rep);

n = 512;
x = [0:n/2,-n/2+1:-1]';
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2+Y.^2);

q = 50;
alist = linspace(1,4,q);

Z = 2*rand(n,n)-1;
S = 1;
for i=1:q
    a = alist(i);
    f = real( ifft2( Z ./ ( S + R.^a ) ) );
    f = f-mean(f(:)); f = f/max(abs(f(:)));
    % draw
    r = 10;
    clf; hold on;
    imagesc(f');
    contour(f',r, 'k');
    colormap(parula(r+1));
    axis off;
    drawnow;
    saveas(gcf, [rep 'Sobolev-' num2str(i) '.png'], 'png');
end
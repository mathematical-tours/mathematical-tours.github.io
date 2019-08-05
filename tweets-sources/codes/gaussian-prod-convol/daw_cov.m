function daw_cov(c, col)


r = 15; % #levelines
m = linspace(0,1,r-1)';

n = 255;
t = linspace(-3,3,n);
[Y,X] = meshgrid(t,t);

ECov  = @(C)C(1,1).*X.^2 + 2*C(1,2).*X.*Y + C(2,2).*Y.^2;
ECovI = @(C)ECov(inv(C));
Gauss = @(C)exp(-ECovI(C)/2);


clf; hold on;
imagesc(t,t,Gauss(c)');
contour(t,t,Gauss(c)',linspace(0,1,r), 'k');
colormap( m*col + (1-m)*[1 1 1] );
caxis([0 1]);
axis equal; box on; set(gca, 'XTick', [], 'YTick', []);

end
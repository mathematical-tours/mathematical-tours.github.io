%%
% conic distance

nx = 200; ny = 150;
x1 = linspace(-1,1,nx);
y1 = linspace(0,1.5,ny);
[Y,X] = meshgrid(y1,x1);

[T,R] = cart2pol(X,Y);

x = 0; y = 1;
[t,r] = cart2pol(x,y);

conic_dist = @(r1,t1,r2,t2)r1.^2 + r2.^2 - 2*r1.*r2.*cos(abs(t1-t2));

d = conic_dist(r,t,R,T);

clf;
imagesc(x1,y1,d');
axis xy;
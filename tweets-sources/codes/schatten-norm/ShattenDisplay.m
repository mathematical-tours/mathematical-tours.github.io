%%
% Display shatten p-norm balls of 
% [a,c/2;c/2,b]

n = 129;
s = 1.4;
t = linspace(-s,s,n);
[A,B,C] = ndgrid(t,t,t);

D = @(A,B,C,s)( (A+B)+s*sqrt((A-B).^2+C.^2) )/2;
S = @(A,B,C,p)abs(D(A,B,C,+1)).^p + abs(D(A,B,C,-1)).^p;

p = 15;

q = 30;
plist = linspace(1,5,q);

for i=1:q
m = (i-1)/(q-1);
col = [m 0 1-m];
p = plist(i);
F = S(A,B,C,p).^(1/p);
clf;
pa = patch( isosurface( t,t,t, F, .6 ) );
isonormals( t,t,t,F,pa );
set(pa, 'FaceColor', col, 'EdgeColor', 'none');
box on; axis off; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
axis equal; axis(s*[-1 1 -1 1 -1 1]);
view(-20,-65);
% zoom(1.1); 
lighting phong;
camlight; drawnow;
end
n = 1000;
[B,A] = meshgrid(1:n,1:n);
F = mod(sqrt(A.^2+B.^2),1);

for i=5:n
clf;
imagesc(F(1:i,1:i)<1e-4);
% imagesc(F(1:i,1:i));
axis equal;
drawnow;
end
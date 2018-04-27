function [I,M] = RenderScene(P,f,S,oc,r)

p = length(f);
n = sqrt(p/6);

% get back the five walls
resh = @(u)reshape(u, [n n size(u,3)]);
extr = @(f,i)resh( f((i-1)*n^2+1:i*n^2,:,:) );
I = {}; M = {};
for i=1:5
    I{i} = extr(f,i);
    M{i} = extr(P,i);
end


% display in 3D
hold on;
for i=1:5
    surf(M{i}(:,:,1), M{i}(:,:,2), M{i}(:,:,3), I{i});
end
% add occluders
for k=1:length(oc)
    a = oc{k};
    [x,y,z] = sphere(30);
    C = x*0; % , x*0, x*0+1; 
    C = cat(3, x*0, x*0, x*0+1);
    h = surf(r(k)*x+a(1),r(k)*y+a(2),r(k)*z+a(2), C);
    alpha(h, .1);
end
% add light source
J = find(S(:)>0); 
PX = P(:,:,1); PY = P(:,:,2); 
ux = [min(PX(J)); max(PX(J))]; uy = [min(PY(J)); max(PY(J))]; 
tx = linspace(ux(1),ux(2), 30); ty = linspace(uy(1),uy(2), 30); 
[y,x] = meshgrid(ty,tx); 
C = cat(3, x*0+1, x*0, x*0); z = x*0+1;
h = surf(x,y,z, C); alpha(h, .1);
%
caxis([min(f(1:5*n*n)) max(f(1:5*n*n))]);
colormap gray(256);
axis ij; view(150,55); % change 150 for rotations
shading interp; axis off;
k = 0;
axis equal;

end
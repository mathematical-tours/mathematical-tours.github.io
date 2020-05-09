addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


dotp = @(x,y)sum(x(:).*y(:));


% [Q,R] = qr(randn(3));
A = diag([1 .94 .93]);
A = diag([1 .94 .8]);
A = diag([1 1 .5]);
[V,D] = eig(A); 
[~,iv] = max(abs(diag(D)));

if 1
V = 0;
while norm(imag(V))<1e-3
A = randn(3);
[V,D] = eig(A); 
end
% norm(A-V*D*inv(V))
[~,I] = sort(diag(abs(imag(D))), 'ascend');
% V = real(V(:,I));
iv = 1;
D = diag( [1 1.06*exp(2i*pi*.01) 1.06*exp(-2i*pi*.01)] );
A = real(V*D*inv(V));
V = real(V);
end

r = 1*1.05;
t = 0.03*pi;

A = [1 0 0; 0 r*cos(t) r*sin(t); 0 -r*sin(t) r*cos(t)];
v = [1 0 0]';

m = 12;
col = distinguishable_colors(m);

u = randn(3,m);
u = u./sqrt(sum(u.^2));
q = 80;
clf; hold on;
plot3(50*[v(1) -v(1)], 50*[v(2) -v(2)], 50*[v(3) -v(3)], 'k', 'LineWidth', 2);
[x,y,z] = sphere(100);
surface(x,y,z, repmat(x*0+1, [1 1 3]));
alpha(.3);
view(35,45);
lighting phong; shading interp; camlight;
for it=1:q
    s = (it-1)/(q-1);
    for k=1:m
        c = col(k,:)*(1-s) + [0 0 0]*s;
        c = col(k,:);
        plot3(u(1,k), u(2,k), u(3,k), '.','color', c, 'MarkerSize', 26);
    end
    axis([-1 1 -1 1 -1 1]); axis square; 
    axis off;
    drawnow;
    u = A*u;
    u = u./sqrt(sum(u.^2));
    mysaveas(it);
end

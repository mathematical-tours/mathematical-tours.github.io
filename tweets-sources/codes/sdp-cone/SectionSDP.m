%%
% Display random section SDP cone.

addpath('../toolbox/');
rep = MkResRep();

d = 20; 

[Q,R] = qr(randn(d*d));
for i=1:4
    A{i} = randn(d); A{i} = A{i}+A{i}'; 
end

n = 301;
t = linspace(-1,1,n)*3;

E = zeros(n);

q = 70;
for it=1:q
    s = (it-1)/(q-1);
    U = s*A{1} + (1-s)*A{2}; U = U/norm(U(:));
    V = s*A{3} + (1-s)*A{4}; V = V/norm(V(:));
    for i=1:n
        for j=1:n
            M = eye(d) + t(i)*U + t(j)*V;
            E(i,j) = min(eig(M));
        end
    end
    I = rescale(-(E>0));
    E = E/max(abs(E(:)));
    % display quantized colormap
    r = 15; % #levellines
    clf; hold on;
    imagesc(t,t,E);
    contour(t,t,E,linspace(-1,1,r), 'k');
    contour(t,t,E,[0 0], 'k', 'LineWidth', 3);
    colormap(parula(r-1));
    caxis([-1 1]);
    axis image; axis off;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end
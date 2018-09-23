% display of levelset with increasing resolution


rep = MkResRep('marching');

name = 'bunny';
[V,Fc] = read_off([name '.off']);

q = 50;
nlist = round(linspace(8,140,q));

for it=1:q
    s = (it-1)/(q-1);
    n = nlist(it);
    %
    t = linspace(-1,1,n);
    [x,y,z] = ndgrid(t,t,t);
    % signed distance function
    V = V - repmat(mean(V,2), [1 size(V,2)]);
    V = V / max(abs(V(:))) * .95;
    A = inpolyhedron(Fc',V',t,t,t);
    D = pointCloud2TDF(V,y,x,z);
    F = -(2*A-1) .* D;
    F = smooth3d(F,4);
    % display
    clf;
    p = patch( isosurface( t,t,t, F, 0 ) );
    isonormals( t,t,t,F,p );
    p.FaceColor = [s 0 1-s];
    axis off;
    axis equal; axis([-1 1 -1 1 -1 1]);
    view(0,90); % view(-150,-40);
    zoom(1.1); lighting phong;
    camlight; drawnow;
    saveas(gcf, [rep 'wireframe-' znum2str(it,2) '.png']);    
    % p.EdgeColor = 'none';
    % saveas(gcf, [rep 'smooth-' znum2str(it,2) '.png']);
end


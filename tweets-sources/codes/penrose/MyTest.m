
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);



% file:///Users/gpeyre/Desktop/penrose/doc/html/PenroseRhombusTiling.html
% Pentaplexity
t = table;
for k = 1:5
    thetad = 72*(k-1);
    t_a = aTriangle(0,cosd(thetad) + 1i*sind(thetad),[]);
    t_ap = apTriangle(0,t_a.Right,[]);
    t = [t ; t_a ; t_ap];
end

t = bTriangle([],-1,1);


kmax = 8;
it=0;
for k = 1:kmax
    for s=1:5
        it = it+1;
        %
        clf;
        showTriangles(t);
        % showDecoratedTiles(t)
        axis equal
        drawnow;
        mysaveas('triangles', it);
        %
        clf;
        showTiles(t);
        % showDecoratedTiles(t)
        axis equal
        drawnow;
        mysaveas('tiles', it);
    end
    
    
    if k<kmax
        t = decomposeTriangles(t);
    end
end
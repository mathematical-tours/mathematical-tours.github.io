
addpath('../toolbox/');
rep = MkResRep('3d');

n = 201;
a = 1.05;
x = linspace(-a,a,n);
[X,Y,Z] = meshgrid(x,x,x);


q = 50; 
plist = linspace(.5,5,q);
vmax = 1;

for i=1:q
    s = (i-1)/(q-1);
    %
    p = plist(i);
    if p~=Inf
        F = ( abs(X).^p + abs(Y).^p+ abs(Z).^p).^(1/p);
    else
        F = max(abs(X),abs(Y));
    end
    %
    % display
    clf;
    p = patch( isosurface( x,x,x,F, vmax ) );
    isonormals( x,x,x,F,p );    
    R = ones(n,n,n);
    isocolors(x,x,x,R*s,R*0,R*(1-s),p);
    p.FaceColor = 'interp';
    p.EdgeColor = 'none';
    box on; axis on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    axis equal; axis([-1 1 -1 1 -1 1]);
    lighting gouraud;
    view(3);    
    camlight; drawnow;  axis off; 
    saveas(gcf, [rep 'balls-' znum2str(i,2) '.png'], 'png');
end


% AutoCrop(rep, ['balls-'], 50); 
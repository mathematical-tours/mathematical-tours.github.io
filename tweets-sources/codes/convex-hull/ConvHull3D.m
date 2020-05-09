% convex hull in 2D

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep  name '-' znum2str(it,3) '.png']);


rand('state', 123);
randn('state', 123);

% roots 
k = 15;
x = rand(k,3);
% speed 
eta = .02;
v = randn(k,3);
v = eta * v./sqrt(sum(v.^2,2));

q = 120;
t = linspace(-100,100,1000);

for it=1:q
    s = (it-1)/(q-1);
    
    ch = convhull(x(:,1),x(:,2),x(:,3));
    I = unique(ch(:));
    clf; hold on;
    plot3(x(:,1), x(:,2), x(:,3), '.', 'color', [s 0 1-s], 'MarkerSize', 10);
    plot3(x(I,1), x(I,2), x(I,3), '.', 'color', [s 0 1-s], 'MarkerSize', 25);
    h = trisurf(ch,x(:,1), x(:,2), x(:,3),'FaceColor',[s 0 1-s],'EdgeColor',[s 0 1-s]);
    h.FaceVertexCData = ones(size(x,1),1) * [s 0 1-s];
    alpha(.8); 
    axis equal; box on;
    view(3);
    shading faceted; camlight;
    e = .15;
    axis([-e 1+e -e 1+e -e 1+e]);
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    drawnow;
    % mysaveas('anim',it);
    
    % advance
    x = x + v;
    % reflexion on boundary
    for d=1:3
        I = find( x(:,d)<0 | x(:,d)>1 );
        v(I,d) = -v(I,d);
    end
end

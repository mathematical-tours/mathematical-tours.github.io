%%
% Display cyclic polytope on the moment curve

q = 50;

for it=4:q
    t = linspace(0,1,it);
    [x,y,z] = deal(t,t.^2,t.^3);
    ch = convhull(x,y,z);
    %
    clf; hold on;
    plot3(x,y,z, '.', 'color', [s 0 1-s], 'MarkerSize', 10);
    trisurf(ch,x,y,z,'FaceColor',[s 0 1-s],'EdgeColor',[s 0 1-s]);
    % alpha(.1);
    axis equal; box on;
    view(3); shading faceted;
    camlight;
    e = .15;
    axis([-e 1+e -e 1+e -e 1+e]);
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    drawnow;
    % mysaveas('anim',it);
end
%%
% 3D Metaballs

addpath('../toolbox/');
rep = MkResRep('3d');

n = 120;


K = 4; % #balls
% position of the centers

name = 'gauss';
name = 'inv';
name = 'invb';
name = 'inv2';
switch name
    case 'inv2'
        phi = @(d)1./(.01+d);
    case 'inv'
        phi = @(d)1./(.01+sqrt(d));
    case 'inva'
        phi = @(d)1./(.01+d.^1.1);
    case 'invb'
        phi = @(d)1./(.01+d.^.3);
    case 'gauss'
        s = .08;
        phi = @(d)exp(-d/(2*s^2));
end
Tp = .05;

rand('seed', 1235);
w = rand(K,1)/2+1/2;  % weights
q = 75; % #movie steps


% box
eta=.15; 
BB = [eta,1-eta,eta,1-eta,eta,1-eta];
% initial centers
c0 = eta + (1-2*eta)*rand(K,3);
% initial speed normalized
s0 = randn(K,3); 
s0 = 4*s0 ./ repmat( sqrt(sum(s0.^2,2)), [1 3] );
% time step
tau = .015;

x = linspace(0,1,n);

c = c0;
s = s0;
g = 0; % gravity constant
for i=1:q
    % generate
    [F,R,T] = GenBall3D(c,w,n, phi, Tp);
    % display
    clf;
    p = patch( isosurface( x,x,x,F, T ) );
    isonormals( x,x,x,F,p );    
    isocolors(x,x,x,R(:,:,:,1),R(:,:,:,2),R(:,:,:,3),p);
    p.FaceColor = 'interp';
    p.EdgeColor = 'none';
    box on; axis on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    axis equal; axis([0 1 0 1 0 1]);
    lighting gouraud;
    view(3);    
    camlight; drawnow;
    % update position and speeds
    c = c + tau*s; 
    % gravity
    s(:,3) = s(:,3) - tau*g;
    % boundary condition
    [c,s] = BoundaryConditions(c,s, BB);
    % save
    saveas(gcf, [rep 'balls' num2str(K) '-' znum2str(i,3) '.png'], 'png');
end


% AutoCrop(rep, ['balls' num2str(K) '-']); %=  to crop the image generated using saveas
% convert balls4-*.png balls4.gif

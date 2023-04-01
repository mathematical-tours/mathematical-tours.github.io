%%
% Planar particules.

addpath('../toolbox/');
rep = MkResRep();

name = 'single';
name = 'two';

save_film = 0;

% #max particles
p = 500;

X = [];
V = [];

tau = .005; % time step for simulation
niter = 2000;
ndisp_tot = 200;
ndisp = round(niter/ndisp_tot);

xmax = 5; % for display
ymax = 3; 

switch name
    case 'single'
        r = 10; % rate of creation
        x0 = {[0;1]}; % initial position
        v0 = {[0;1]}; % average initial speed
        sigma = .25; % dispertion
    case 'two'
        r = 12*4; % rate of creation
        x0 = {[-2;0], [+2;0]}; % initial position
        v0 = {[.5;.9]*2, [-.5;.9]*2}; % average initial speed
        sigma = .15; % dispertion
end

g = -1*[0;1]; % gravity
eta = .8; % damping of rebound
kappa = .3; % radius of the balls

Col = distinguishable_colors(p);
save_film = 1;

ColM = { {[1 0 0], [1 1 0]}, {[0 1 1] [0 0 1]} };
isvg = 1;
for i=1:niter
    % add new particles
    if mod(i,r)==1
        % select with source
        m = ceil(rand*length(v0));
        % advect
        X(:,end+1) = x0{m}; 
        V(:,end+1) = v0{m} + randn(2,1)*sigma;
        n = size(X,2);
        %
        h = rand;
        Col(n,:) = h*ColM{m}{1} + (1-h)*ColM{m}{2};
    end
    % evolve
    X = X + tau*V;
    V = V + tau*g*ones(1,size(X,2));
    % kill extra particules
    k = max(size(X,2)-p,0);
    X(:,1:k) = [];  V(:,1:k) = [];
    Col = [Col(k+1:end,:); Col(1:k,:)]; % shift colors
    n = size(X,2);
    % colision with ground
    I = find(X(2,:)<0);
    X(2,I) = 0; V(2,I) = - eta * V(2,I);
    % self collision
    D = distmat(X,X);
    [B,A] = meshgrid(1:n,1:n);
    I = find(D<kappa & B>A);
    A = A(I); B = B(I);
    % treat each collision 
    for k=1:length(B)
        a = A(k); b = B(k);
        V(:,[a b]) = V(:,[b a]);
    end
    % display
    s = ones(n,1)*100; % size
    clf; hold on; 
    scatter( X(1,:), X(2,:), s, Col(1:n,:), 'filled' );
    axis equal; axis off;
    plot([-xmax xmax], [0 0]-kappa/2, 'k', 'LineWidth', 3);
    axis([-xmax,xmax,-kappa,ymax]); drawnow;
    if mod(i,ndisp)==1 && save_film
        saveas(gcf, [rep name '-' znum2str(isvg,3) '.png']);
        isvg = isvg+1;
    end
end



% AutoCrop(rep, [name '-']); 
% convert balls4-*.png balls4.gif

%%
% Verlet integrations. 

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

xmax = 2;
g = [0 -.5]';

x0 = [.1 .5]';
v0 = [1 1]'*.3;

k = 15; % #particles
x0 = rand(2,k); 
x0(1,:) = x0(1,:)*xmax;
%
v0 = randn(2,k)*.2; 

radius = .12;
radius = .08;

tau = .01/2;
sw = 30;

tau = .01/2;
sw = round(10);



Col = distinguishable_colors(k);

X = x0; v = v0;
niter = 2000;
it1 = 0;
for it=1:niter
    if mod(it,sw)==1
        clf; hold on;
        for j=1:k
            mem = 200;
            plot(squeeze(X(1,j,max(1,end-mem):end)), squeeze(X(2,j,max(1,end-mem):end)), 'Color', Col(j,:), 'LineWidth', 2);
            plot(X(1,j,end), X(2,j,end), 'k.', 'MarkerSize', 580*radius, 'Color', Col(j,:));
        end
        set(gca, 'XTick', [], 'YTick', []);
        axis equal;
        axis([0 xmax 0 1]); box on; %
        drawnow;
        it1 = it1+1;
        mysaveas(it1);
    end
    % verlet update
    if 1
        v = v + tau*g;
        x = X(:,:,end) + tau*v;
    else
        [v,x] = deal(v + tau*g, X(:,:,end) + tau*v);
    end
    X(:,:,end+1) = x;
    % bounding
    I = find(X(2,:,end)<0);
    X(2,I,end) = 0;
	v(2,I) = -v(2,I);
    %    
    I = find(X(1,:,end)<0);
    X(1,I,end) = 0;
	v(1,I) = -v(1,I);
    %    
    I = find(X(1,:,end)>xmax);
    X(1,I,end) = xmax;
	v(1,I) = -v(1,I);
    % pairwise collision
    D = distmat(X(:,:,end), X(:,:,end));
    [a,b] = find(sparse(D<radius));
    [a,b] = deal(a(a<b), b(a<b));
    [v(:,a),v(:,b)] = deal(v(:,b),v(:,a));
end


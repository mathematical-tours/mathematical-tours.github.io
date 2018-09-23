%%
% Test for 2D meshing.

addpath('../toolbox/');
rep = MkResRep();

% draw polygon


if 0
    % click selection for outside/inside
    it = 0;
    node = [];
    clf; hold on;
    K = [];
    for k=1:1
        K(end+1) = 0;
        while true
            if size(node,1)>1
                if length(K)==1
                    plot(node(:,1), node(:,2), 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
                else
                    plot(node([1:K(1) 1],1), node([1:K(1) 1],2), 'b.-', 'LineWidth', 2, 'MarkerSize', 20);
                    plot(node(K(1)+1:end,1), node(K(1)+1:end,2), 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
                end
            end
            axis([0 1 0 1]);
            [a,b,button] = ginput(1);
            plot(a,b, '.', 'MarkerSize', 15);
            if button==3
                break;
            end
            node(end+1,:) = [a b]';
            K(end) = K(end)+1;
        end
    end
    
    n = size(node,1);
    edge1 = [1:K(1); [2:K(1) 1]]';
    edge2 = [K(1)+1:n; [K(1)+2:n K(1)+1]]';
    edge = [edge1]; % ;edge2];
else
    
%% load from a binary image
f = imread('france.jpg');
f = sum(f,3); f = double((f/max(f(:)))>.5);
[C,h] = contour(f,[.5,.5]);
m = C(2,1);  c = C(:,2:m+1); c = c(1,:)'+1i*c(2,:)';
c = (c-1)/(max(size(f))-1);
node = [real(c), imag(c)];

% meshfile = './poly-data/lake.msh';
% [node,edge] = triread( meshfile );
end

% curve re-sampling
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );



node0 = node;

q = 20;
pblist = round(linspace(50,1500,q));
for it=1:length(hlist)
    
    % re-sample
    p = pblist(it);
    c = node0(:,1) + 1i*node0(:,2);
    c = resample(c,p);
    node = [real(c), imag(c)];
    n = size(node,1);
    edge = [1:n; [2:n 1]]';
    
    
    
    % generate
    opt.disp = Inf;
    opt.rho2 = 1.1;
    [vert,etri, tria,tnum] = refine2(node,edge,[],opt); % hfun
    % polish
    [vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum,opt);
    % draw
    clf;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    axis([0 1 0 1]);
    axis ij;
    drawnow;
    saveas(gcf, [rep 'mesh2d-' znum2str(it,2) '.png']);
end
% AutoCrop(rep, 'mesh2d-');
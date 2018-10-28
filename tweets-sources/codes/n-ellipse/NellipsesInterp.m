%%
% Draw levelsets for n-ellipses.


if not(exist('test'))
    test=1;
end

addpath('../toolbox/');
rep = MkResRep(num2str(test));

% click and play
clf; hold on;
Q = {};
for k=1:2
    Q{k} = [];
    while true
        axis([-1 1 -1 1]);
        [a,b,button] = ginput(1);
        if button==3
            break;
        end
        plot(a,b, '.', 'MarkerSize', 15);
        Q{k}(:,end+1) = [a;b];
    end
end

N = 256;
x = linspace(-1,1,N);
[Y,X] = meshgrid(x,x);

% #points
P = min(size(Q{1},2),size(Q{2},2));
Q{1} = Q{1}(:,1:P);
Q{2} = Q{2}(:,1:P);

q = 50;
for it=1:q
    
    t = (it-1)/(q-1);
    Qi = Q{1}*(1-t) + Q{2}*t;
    
    D = zeros(N);
    for i=1:P
        D = D + sqrt( (X-Qi(1,i)).^2 + (Y-Qi(2,i)).^2 );
    end
    
    D1 = sqrt(rescale(D));
    r = 12;
    clf; hold on;
    imagesc(x,x,D1);
    contour(x,x,D1, linspace(0,1,r), 'k', 'LineWidth', 2);
    plot(Qi(2,:), Qi(1,:), 'r.', 'MarkerSize', 30);
    axis image; axis off;
    colormap( parula(r-1) );
    caxis([0 1]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    
end

% AutoCrop(rep, 'anim-');

test = test+1;

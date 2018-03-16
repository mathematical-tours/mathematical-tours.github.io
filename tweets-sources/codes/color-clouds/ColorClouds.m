
if not(exist('name'))
name = 'picasso';
end

rep = ['../results/color-clouds/' name '/'];
[~,~] = mkdir(rep);
addpath('../toolbox/');

n = 200;
f = load_image(name, n);
f = rescale(f);

% planar display
m = 50;

for it=1:m
    t = (it-1)/(m-1);
    t = 1e-2 + t^2;
    
    if 0
    ms = 25;
    clf; hold on;
    for i=1:n
        for j=1:n
            z = 1-(i-1)/(n-1); x = (j-1)/(n-1); y = 1/2;
            c = f(i,j,:); c = c(:);
            h = plot3(...
                (1-t)*x + t*c(1), ...
                (1-t)*y + t*c(2), ...
                (1-t)*z + t*c(3), ...
                '.', 'MarkerSize', ms, 'Color', c);
        end
    end
    end
    
    % base color
    c = reshape(f, n*n,3);
    % base position
    r = linspace(0,1,n);
    [X,Z] = meshgrid(r,1-r); Y = zeros(n,n)+1/2; 
    %
    s = ones(n*n,1)*10; % size
    clf;
    scatter3(...
        t*c(:,1) + (1-t)*X(:), ...
        t*c(:,2) + (1-t)*Y(:), ...
        t*c(:,3) + (1-t)*Z(:),s, c, 'filled' );
    %    
    axis ij; view(3); drawnow;
    axis equal; axis([0 1 0 1 0 1]);
    box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    saveas(gcf, [rep name '-' num2str(it) '.png']);
end
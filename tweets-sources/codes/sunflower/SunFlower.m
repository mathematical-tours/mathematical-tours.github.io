%%
% Display a sunflower spiral with voronoi cells
% https://www.irishtimes.com/news/science/the-efficient-use-of-space-behind-ibec-s-sunflower-style-logo-1.1820736
% https://thatsmaths.com/2014/06/05/sunflowers-and-fibonacci-models-of-efficiency/
% Vogel, H (1979). "A better way to construct the sunflower head". Mathematical Biosciences. 44 (44): 179?189


addpath('../toolbox/');
rep = MkResRep();

theta = pi/12;

% golden angle
theta = pi*137.5/180;

theta = pi/10;

n = 50; %# points
r = 1:n-1;


q = 90;
t_list = linspace(.5,.8,q)*pi;

t_list = 137.5/180 * linspace(.99,1.01,q)*pi;

for it=1:q
    theta = t_list(it);
    s = (it-1)/(q-1);
    % generate
    Z = .95*sqrt(r/n).*exp( 1i*r*theta );
    % voronoi
    [VX,VY] = voronoi(real(Z), imag(Z));
    V = VX + 1i*VY;   
    
    % Display Delaunay
    clf; hold on;
    plot(V, '-', 'Color', [s 0 1-s], 'LineWidth', 2);
    plot(Z, '.', 'Color', [s 0 1-s]*.8,  'MarkerSize', 25);
    axis equal;axis equal; axis([-1 1 -1 1]);
    axis on; box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'anim');
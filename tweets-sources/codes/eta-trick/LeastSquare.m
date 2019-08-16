% Do robust least squares.

% #points
n = 40;
% data, well spread
t = cumsum( .5 + rand(n,1) );
t = rescale(t, .05, .95);

%
y = rand(n,1);
%
y = (t-1/2).^3 + .03 * randn(n,1);

% degree for fitting
d = 6;

dmax = 50;
for d=0:dmax
    
    % features
    x = t .^( 0:d );
    %
    % for display
    T = linspace(0,1,512)';
    X = T .^( 0:d );
    
    % L2
    % min |y-X*w|
    w = pinv(x)*y;
    %
    clf; hold on;
    plot(t, y, '.k', 'MarkerSize', 15);
    plot(T, X*w, 'color', [d/dmax 0 1-d/dmax], 'LineWidth', 2);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'XTick', [], 'YTick', []);
    axis([0 1 -1/8 1/8]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(d+1,2) '.png']);
end
%%
% Plot of Weierstrass function.



addpath('../toolbox/');
rep = MkResRep();

n = 512*10; 
K = 200*2;


% a>1
% ab>1
q = 50;
t = (0:q-1)/q;


a_list = t*0 + .6;
b_list = .2*(1-t) + 2*t; 
x = linspace(0,2*pi,n);
vmin = -2; vmax = 4;


a_list = .2*(1-t) + .8*t;
b = 2;
b_list = b*(1-t) + b*t; 
x = linspace(0,1,n);
vmin = -2; vmax = 2;


[K,X] = meshgrid(0:K-1,x);

W = @(a,b)sum( (a.^K) .* cos( (b.^K) .* pi .* X), 2 );




    
for i=1:q
    a = a_list(i);
    b = b_list(i);
    t = (i-1)/(q-1);
    y = W(a,b);
    y = rescale(y,.1, .9);
    clf;  hold on;    
    area(x, y, 0,  'FaceColor', .5+.5*[1-t 0 t]);
    plot(x, y, 'color', [1-t 0 t], 'LineWidth', 2);
    axis tight; axis([min(x) max(x) 0 1]); % vmin vmax
    box on;
    % axis off;
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'evol-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, 'evol-');
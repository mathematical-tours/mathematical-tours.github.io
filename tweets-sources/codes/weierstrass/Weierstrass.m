%%
% Plot of Weierstrass function.


rep = '../results/weierstrass/';
[~,~] = mkdir(rep);

n = 512*10; 
K = 200*2;

x = linspace(0,1,n);
[K,X] = meshgrid(0:K-1,x);

W = @(a,b)sum( (a.^K) .* cos( (b.^K) .* pi .* X), 2 );

% a>1
% ab>1
q = 50;
t = (0:q-1)/q;

a_list = .2*(1-t) + .8*t;
b_list = 2*(1-t) + 2*t; 

vmin = -2; vmax = 2;
    
for i=1:q
    a = a_list(i);
    b = b_list(i);
    t = (i-1)/(q-1);
    clf;  hold on;
    area(x, W(a,b), vmin,  'FaceColor', .5+.5*[1-t 0 t]);
    plot(x, W(a,b), 'color', [1-t 0 t], 'LineWidth', 2);
    axis tight; axis([min(x) max(x) vmin vmax]);
    axis off;
    drawnow;
    saveas(gcf, [rep 'evol-' num2str(i) '.png']);
end
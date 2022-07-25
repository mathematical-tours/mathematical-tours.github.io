addpath('../toolbox/');
global rep
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 1400;
n = 1000/2;
noise = 0;
global param_t;

name = 'swiss';
name = 'helix';
[X, labels, param_t] = generate_data(name, n, noise);
X(:,2) = X(:,2)*.7;

clf;
scatter3(X(:,1),X(:,2),X(:,3), 20, param_t(:,1), 'filled');
box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
axis equal; 
saveas(gcf, [rep 'cloud.png']);

Y0 = randn(n,2)*2;



q = 150;
MaxIter = 15000;
global disp_list
global cnt 
disp_list = unique( round( 1 + (MaxIter-1)*linspace(0,1,q).^4 ) );
cnt = 1;

perp = 15;

% 'InitialY', Y0, 
opts = statset('OutputFcn',@TSNE_Callback, 'MaxIter', MaxIter);
Y = tsne(X,'Algorithm', 'exact', 'Standardize',true, 'Perplexity',perp, 'NumPrint', 1, 'Options', opts);

return

q = 80;
perp_list = 1 + 100 * linspace(0,1,q).^3;

for it=1:q
    perp = perp_list(it);
    rng('default');
    % 
    % 'barneshut'
    Y = tsne(X,'Algorithm', 'exact', 'Standardize',true, 'Perplexity',perp, 'InitialY', Y0);
    % Plot the result
    clf;
    scatter(Y(:,1),Y(:,2), 20, t(:,1), 'filled');
    axis equal;
    axis off;
    drawnow;
    mysaveas(it);
end
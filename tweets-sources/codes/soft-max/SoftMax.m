%%
% Plot of soft-max

rep = '../results/soft-max/';
[~,~] = mkdir(rep);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);
rescale = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));


% stabilized log-sum-exp and soft max
max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);
LSE = @(S)log( sum(exp(S), 2) );
LSE = @(S)LSE( S-max2(S) ) + max(S,[],2);
SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);
SM = @(S)SM(S-max2(S));

q = 200;
t = linspace(0,1,q);
[Y,X] = meshgrid(t,t);
U = [X(:), Y(:)];

n = 12;
a = randn(1,n);
eps_list = [1 .75 .5 .25 .1 .01];
for i=1:length(eps_list)
    epsilon = eps_list(i);
    c = (i-1)/(length(eps_list)-1);
    clf; hold on;
    bar(1:n,SM(a/epsilon), 'FaceColor', [c 0 1-c]);
    plot(rescale(a), 'k.--', 'MarkerSize', 25, 'LineWidth', 1);
    axis([.5 n+.5 0 1]);
    box on; set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/3); box on;
    saveas(gcf, [rep 'sm-' num2str(epsilon) '.eps'], 'epsc');
end

eps_list = [1 .5 .1 .01 .001];
for i=1:length(eps_list)
    epsilon = eps_list(i);
    R = reshape(LSE(U/epsilon)*epsilon,[q q]);
    % display
    r = 15; % #levellines
    clf; hold on;
    imagesc(t,t,R');
    contour(t,t,R',linspace(0,1,r), 'k');
    colormap(parula(r-1));
    caxis([0 1]);
    axis image; axis off;
    saveas(gcf, [rep 'lse-' num2str(epsilon) '.png']);
end



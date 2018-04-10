
addpath('../toolbox/');
rep = MkResRep();



% save for sampling 
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));
k = 50; 
U = clamp(randn(k),-3,3);
imwrite( rescale(Upsc(U,5)), [rep 'matrix-gaussian.png'] );
U = randn(k)>0;
imwrite( rescale(Upsc(U,5)), [rep 'matrix-bernoulli.png'] );

ms = 20;
n = 1000;
t = linspace(0,1,1000);

myrand = @(n)randn(n);
myrand = @(n)(2*(randn(n)>0)-1);


%%
% randn
nlist = [100 500 1000 2000];
for k=1:length(nlist)
    n = nlist(k);
    A = myrand(n)/sqrt(n);
    clf; hold on;
    plot( eig(A), '.', 'MarkerSize', ms  );
    plot( cos(2*pi*t), sin(2*pi*t), 'r', 'LineWidth', 2 )
    m = 1.05;
    axis([-m m -m m]); box on;
    axis tight;
    axis equal; box on;
    set(gca, 'FontSize', 20, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'circlelaw-' num2str(n) '.eps'], 'epsc');
end
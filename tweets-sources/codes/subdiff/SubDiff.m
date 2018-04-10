%%
% Just some simple plots on the subdifferential.

addpath('../toolbox');
rep = '../results/subdiff/';
[~,~] = mkdir(rep);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

n = 1025;
t = linspace(-1,1,n);
h = (max(t)-min(t))/(n-1);
mdiff = @(f)( f(2:end) - f(1:end-1) ) / h;
mmean = @(f)( f(2:end) + f(1:end-1) )/2;

F = {.3*abs(t-.5) + .7*abs(t+.5) + .5*t.^2, abs(t)};
names = {'piecereg' 'abs'};

for k=1:length(F)
    name = names{k};
    f = F{k}; 
    f = rescale(f, -1, 1);
    clf;
    plot(t, f, 'b-', 'LineWidth', 2);
    axis normal; axis([-1 1 min(f)-.1 max(f)+.1]);
    set(gca, 'FontSize', 20);
    set(gca, 'XTick', [-1 0 1], 'YTick',  -3:1:3);
    SetAR(1/2);
    saveas(gcf, [rep name '.eps'], 'epsc');
    clf;
    df = mdiff(f);
    plot(mmean(t), df, 'r-', 'LineWidth', 2);
    axis normal; axis([-1 1 min(df)-.1 max(df)+.1]);
    set(gca, 'FontSize', 20);
    set(gca, 'XTick', [-1 0 1], 'YTick', -3:1:3);
    SetAR(1/2);
    saveas(gcf, [rep name '-diff.eps'], 'epsc');
end

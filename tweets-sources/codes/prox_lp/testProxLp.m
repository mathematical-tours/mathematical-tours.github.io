%%
% Plot prox of l^p norm


rep = ['results/'];
[~,~] = mkdir(rep);

y = linspace(-3,3,1025)';

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

tau = 1;
p=1;



pmode = 2;

switch pmode
    case 1
        plist = 1:1/5:2;
    case 2
        plist = 0:1/5:1; plist(1)=1e-3;
end


clf; hold on;
for i=1:length(plist)
    c = (i-1)/(length(plist)-1);
    switch pmode
        case 1
            col = [c 0 1-c];
        case 2
            col = [0 1-c c];
    end
    plot(y, ProxLP(y,plist(i),tau), 'Linewidth', 2, 'Color', col);
end
axis tight;
box on; set(gca, 'FontSize', 20);
SetAR(2/3);
saveas(gcf, [rep 'Prox-' num2str(pmode) '.eps'], 'epsc');

clf; hold on;
for i=1:length(plist)
    c = (i-1)/(length(plist)-1);
    switch pmode
        case 1
            col = [c 0 1-c];
        case 2
            col = [0 1-c c];
    end
    plot(y, abs(y).^plist(i), 'Linewidth', 2, 'Color', col);
end
axis tight;
box on; set(gca, 'FontSize', 20);
SetAR(2/3);
saveas(gcf, [rep 'Lp-' num2str(pmode) '.eps'], 'epsc');
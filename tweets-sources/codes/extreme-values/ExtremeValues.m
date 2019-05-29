%%
% Plot of "Generalized extreme value distribution"

addpath('../toolbox/');
rep = MkResRep();

n = 2048;
A = 4;
s = linspace(-A,A,n);

q0 = 12;
xilist0 = linspace(-.7,.9,q0);
q = 50; % for animation
xilist = linspace(-.7,.9,q);

for j=1:q
    t = (j-1)/(q-1);
    % background display
    clf; hold on;
    for i=1:q0
        t0 = (i-1)/(q0-1);
        plot(s(1:n-1), ExtVal(s,xilist0(i)), 'LineWidth', 1, 'Color', .6*[1 1 1] + .4*[t0 0 1-t0]);
    end
    % display
    plot(s(1:n-1), ExtVal(s,xilist(j)), 'LineWidth', 3, 'Color', [t 0 1-t]);
    set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20);
    axis([-A A 0 .6]); box on;
    drawnow;
    saveas(gcf, [rep 'extreme-' znum2str(j,2) '.png']);
end


 % AutoCrop(rep, 'extreme');

%axis tight;
% set(gca, 'XTick', -A:2:A, 'YTick', 0:.2:1, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20);
% axis([-A A 0 .6]); box on;




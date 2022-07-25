addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 1024;

t = linspace(-1,1,n) * 8;
fscale = 3;
gabor = @(mx,mf,s)real( exp(-(t-mx).^2 / (2*s.^2)) .* exp(1i*mf* fscale * t) );

S = [1/1.3,1.3];
mx_max = 4;
mf_max = 4;

% 1/S

q = 100;

for it=1:q
    u = (it-1)/q;
    v = sin(pi*u).^2;
    mx = cos(2*pi*u)*mx_max;
    mf = sin(2*pi*u)*mf_max;
    s  = (1-v)*S(1) + v*S(2);  
    col = [v, 0, 1-v]; 

    clf;
%    subplot(2,1,1);
    plot(t, gabor(mx,mf,s), 'color', col, 'LineWidth', 2);
    axis([min(t) max(t) -1.05 1.05]);
    set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    mysaveas('space', it);
%    subplot(2,1,2);

    clf;
    plot(t, gabor(mf,mx,1/s), 'color', col, 'LineWidth', 2);
    axis([min(t) max(t) -1.05 1.05]);
    set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    mysaveas('freq', it);
    %drawnow;

    clf;
    sc = 2;
    rectangle('Position',[mx-sc*s, mf-sc/s, 2*sc*s, 2*sc/s ], 'FaceColor', col, 'EdgeColor', col);
    axis equal; 
    axis([-1 1 -1 1]*max(t));
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('box', it);

    
end
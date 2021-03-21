addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

x = linspace(-1,1,1024)*4;
p = 50;
nulist = linspace(.01,5,p);

nulist = .05 + 8*linspace(0,1,p).^(3);

ftype = 'fourier';
ftype = 'kernel';

for i=1:p
    clf; hold on;
    for j=1:i
        s = (j-1)/(p-1);
        nu = nulist(j);
        switch ftype
            case 'kernel'
                y = abs(x).^nu .* besselk(nu,sqrt(2*nu)*abs(x));
            case 'fourier'
                y = (nu + abs(x).^2).^(-(nu+1/2));
                
        end        
        y = y/max(y);
        if j==i
            col = [s 0 1-s];
        else
            col = .3*[s 0 1-s]+.7;
        end
        plot(x,y, 'color', col, 'LineWidth', 1+2*(j==i));
    end
    axis([min(x) max(x) 0 1.02]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(ftype,i)
end


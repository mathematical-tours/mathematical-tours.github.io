addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 2048*8+1;
x = 3*linspace(-1,1,n)';

omega = [0:n/2, -n/2+1:-1]' / n;
omega = [0:(n-1)/2, -(n-1)/2:-1]' / n;

lapl = @(x,p) real(ifft( abs(omega).^p .* fft(x) ) );
lapl = @(x,p) real(ifft( abs(sin(omega*pi)).^p .* fft(x) ) );

s = .06;
f = exp(-x.^2/(2*s^2));

clf;
plot(x,f, 'k', 'LineWidth',2);
u = .3;
axis([-u u -.05 1.05]);
box on; set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
drawnow;
saveas(gcf, [rep 'function.png']);

q = 100;
plist = linspace(0,6,q);

for it=1:q
    clf; hold on;
    for jt=1:it
        s = (jt-1)/(q-1);
        col = [s 0 1-s];
        p = plist(jt);
        g = lapl(f,p);
        g = g/max(g);
        plot(x,g, 'color', .5 + .5*col, 'LineWidth', 1);
    end    
    plot(x,g, 'color', col, 'LineWidth', 3);
    axis tight; 
    axis([-u u -.75 1.02]);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end


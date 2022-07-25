addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 512*2;
x = linspace(0,1,n);

omega = [0:n/2, -n/2+1:-1]' * (2/n);


HeatSol = @(f,p,t) real(ifft( exp(-t*abs(omega).^p) .* fft(f) ) );

f0 = zeros(n,1);
f0(round([.3 .55 .7]*end)) = [.3 .7 -1];
f0 = cumsum(f0);






p = .2;
tmax = 8;
texp = 1;


p = 2;
tmax = 300000;
texp = 1;



p = 3.5;
tmax = 200*3000000;
texp = 2.5;

q = 100;
tlist = tmax * linspace(0,1,q).^texp ;
for it=1:q     
    clf; hold on;
    for jt=1:it
        s = (jt-1)/(q-1);
        col = [s 0 1-s];
        g = HeatSol( f0,p,tlist(jt) );
        plot(x,g, 'color', .5 + .5*col, 'LineWidth', 1);
    end    
    plot(x,g, 'color', col, 'LineWidth', 3);
    axis([0 1 -.1 1.1]);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it); 
end
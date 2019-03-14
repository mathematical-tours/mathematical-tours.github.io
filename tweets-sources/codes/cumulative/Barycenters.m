addpath('../toolbox/');
rep = MkResRep('barycenters');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
SetTickOff = @()set(gca, 'XTick',[], 'YTick',[]);
SetTickOn = @()set(gca, 'XTick',[0 1/2 1], 'YTick',[0 1/2 1]);

N = 256;
t = linspace(0,1,N)';
gauss = @(m,s)exp(-(t-m).^2/(2*s).^2);
normalize = @(x)x/sum(x(:));

s = .6;
A = @(t)gauss(.3*s+t,.05*s) + .5*gauss(.6*s+t,.15*s);
B = @(t).5*gauss(1-s*.2-t,.04*s) + .8*gauss(1-s*.5-t,.05*s) + .5*gauss(1-s*.8-t,.04*s);
vmin = .025;
A = @(t)normalize(vmin + A(t));
B = @(t)normalize(vmin + B(t));

a = A(0);
b = B(0);

plot(t,[a b]);

clf;
area(t, a, 'FaceColor', 'r', 'EdgeColor', 'r');
axis tight; SetAR(1/2); SetTickOff();
saveas(gca, [rep 'input-1.png']);

clf;
area(t, b, 'FaceColor', 'b', 'EdgeColor', 'b');
axis tight; SetAR(1/2); SetTickOff();
saveas(gca, [rep 'input-2.png']);


% cumulative
ca = cumsum(a);
cb = cumsum(b);
% inverse cumulatives
ica = interp1(ca, t, t, 'spline');
icb = interp1(cb, t, t, 'spline');

q = 50;

% Barycenter of inverse cumulant`
for it=1:q
    r = (it-1)/(q-1);
    icm = (1-r)*ica + r*icb;
    cm = interp1(icm, t, t, 'spline');
    m = diff([0;cm]);
        
    col = [1-r 0 r];
    
    clf; hold on;
    plot(t, ica, 'color', .2*[1 0 0] + .8*[1 1 1], 'LineWidth', 2);
    plot(t, icb, 'color', .2*[0 0 1] + .8*[1 1 1], 'LineWidth', 2);
    plot(t, icm, 'color', col, 'LineWidth', 2);
    axis tight; SetAR(1); SetTickOff(); box on;
    saveas(gca, [rep 'interp-icumul-' znum2str(it,2) '.png']);

    clf; hold on;
    plot(t, ca, 'color', .2*[1 0 0] + .8*[1 1 1], 'LineWidth', 2);
    plot(t, cb, 'color', .2*[0 0 1] + .8*[1 1 1], 'LineWidth', 2);
    plot(t, cm, 'color', col, 'LineWidth', 2);
    axis tight; SetAR(1); SetTickOff(); box on;
    saveas(gca, [rep 'interp-cumul-' znum2str(it,2) '.png']);
    
    clf;
    area(t, m, 'FaceColor', col, 'EdgeColor', col);
    axis tight; SetAR(1/2); SetTickOff();
    drawnow;
    saveas(gca, [rep 'bary-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'transports-')


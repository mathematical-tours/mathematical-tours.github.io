% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 512; 
t = linspace(0,1,n)';

f = (.1+t) .* cos(8*pi*t.^3)+t.^2;

clf
plot(t,f);

dmax = 26;

k = 1;
k = 8;
deg_list = 0:k:k*dmax;

sub = 6;
q = sub*dmax;
it = 0;
for i=1:length(deg_list)-1
    % approximation using i and i+1
    P0 = t .^ deg_list(1:i);
    P1 = t .^ deg_list(1:i+1);
    f0 = P0*pinv(P0)*f;
    f1 = P1*pinv(P1)*f;
    for j=1:sub
        it = it+1;
        s = (j-1)/sub;
        fs = (1-s)*f0 + s*f1;
        col = (it-1)/(q-1);
        clf; hold on;
        plot(t,f, 'k', 'LineWidth', 2);
        plot(t,fs, 'color', [col,0,1-col], 'LineWidth', 2);
        axis([0, 1, min(f)-.2, max(f)+.2]);
        set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
        box on;
        drawnow;
        mysaveas(it);
    end
end
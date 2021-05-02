%%
% Trivial NN.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

phi = @(t)max(t,0);
f = @(x,a,b,c)c(1)*phi(a(1)*x+b(1)) + c(2)*phi(a(2)*x+b(2));

x = linspace(-1,1,200) * 1.3;
c = [1 -1];

u = .3; v = .6;
a = [1 1] * 1/(v-u);
b1 = -u*a(1);
b = [b1 1-b1];

q = 120;
bmax = 4;
amax = 4;
cmax = 2;
for it=1:q
    t = (it-1)/(q-1);
    u = cos(2*pi*t); v = sin(2*pi*t);
    %
    a = [1 1] * 1/(v-u);
    b1 = -u*a(1);
    b = [b1 b1-1];  
    %
    clf; hold on;
    if 1
        plot(x, f(x,a,b,c), 'LineWidth', 3 );
        axis([min(x) max(x) -.15 1.15]);
        box on;
        set(gca, 'XTick', [], 'YTick', []);
    else
        clf; hold on;
        plot_neuron([0 1], [1 1], b(2)/bmax);
        plot_neuron([0 1], [1 0], b(1)/bmax);
        %
        plot_neuron([0 1], [0 0], a(1)/amax);
        plot_neuron([0 1], [0 1], a(2)/amax);
        %
        plot_neuron([1 2], [0 .5], c(1)/cmax);
        plot_neuron([1 2], [1 .5], c(2)/cmax);
        axis equal;
        box on;
        axis([-.2 2.2 -.2 1.2]);
        set(gca, 'XTick', [], 'YTick', []);
    end
    drawnow;
    mysaveas(it);
end

return;



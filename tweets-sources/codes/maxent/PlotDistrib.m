%%
% Just display some maxent famillies


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = ['../results/maxent/'];
[~,~] = mkdir(rep);

q = 15; % # distrib

n = 1024;
t = linspace(-1,1,1024);
r = linspace(0,1,1024);
Z = @(x)x/sum(x);

% uniforms
alist = linspace(.1,.9,q);
clf; hold on;
Pmax = 0;
for i=1:q
    c = (i-1)/(q-1);
    P = Z(abs(t)<=alist(i));
    Pmax = max(Pmax,max(P));
    plot(t, P,'color', [c 0 1-c], 'LineWidth', 2 );
end
axis([-1 1 -Pmax*.02, Pmax*1.02]);
box on; set(gca, 'XTick', [-1 1], 'YTick', []); set(gca, 'FontSize', 15);
SetAR(1/2);
saveas(gcf, [rep 'uniform.eps'], 'epsc');

% gaussians
mlist = linspace(-.7,.9,q);
slist = linspace(.1,.03,q);
clf; hold on;
Pmax = 0;
for i=1:q
    c = (i-1)/(q-1);
    P = Z( exp(-(t-mlist(i)).^2 / (2*slist(i).^2) ) );
    Pmax = max(Pmax,max(P));
    plot(t, P,'color', [c 0 1-c], 'LineWidth', 2 );
end
axis([-1 1 -Pmax*.02, Pmax*1.02]);
box on; set(gca, 'XTick', [-1 1], 'YTick', []); set(gca, 'FontSize', 15);
SetAR(1/2);
saveas(gcf, [rep 'gaussian.eps'], 'epsc');


% lognormal
mlist = linspace(-3,-.1,q);
slist = linspace(.5,.06,q);
clf; hold on;
Pmax = 0;
for i=1:q
    c = (i-1)/(q-1);
    P = Z( exp(-(log(r)-mlist(i)).^2 / (2*slist(i).^2) ) );
    Pmax = max(Pmax,max(P));
    plot(r, P,'color', [c 0 1-c], 'LineWidth', 2 );
end
axis([0 1 -Pmax*.02, Pmax*1.02]);
box on; set(gca, 'XTick', [0 1], 'YTick', []); set(gca, 'FontSize', 15);
SetAR(1/2);
saveas(gcf, [rep 'lognormal.eps'], 'epsc');


% gamma
tlist = linspace(.1,.2,q);
klist = linspace(1.1,3,q);
clf; hold on;
Pmax = 0;
for i=1:q
    c = (i-1)/(q-1);
    P = Z( r.^(klist(i)-1) .* exp(-r/tlist(i)) );
    Pmax = max(Pmax,max(P));
    plot(r, P,'color', [c 0 1-c], 'LineWidth', 2 );
end
axis([0 1 -Pmax*.02, Pmax*1.02]);
box on; set(gca, 'XTick', [0 1], 'YTick', []); set(gca, 'FontSize', 15);
SetAR(1/2);
saveas(gcf, [rep 'gamma.eps'], 'epsc');

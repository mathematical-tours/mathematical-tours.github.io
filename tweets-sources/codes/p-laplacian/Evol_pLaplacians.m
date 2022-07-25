addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 512/2;
x = linspace(0,1,n);

delta = (max(x)-min(x))/n;
dx  = @(x)1/delta*(x([2:end 1])-x);
dxS = @(x)1/delta*(x-x([end 1:end-1]));
epsilon = 1e-10;
abs1 = @(x)sqrt(x.^2 + epsilon^2);
lapl = @(x,p)dxS( dx(x) .* abs1(dx(x)).^(p-2) );


f0 = zeros(n,1);
f0(round([.3 .55 .7]*end)) = [.3 .7 -1];
f0 = cumsum(f0);





p = 1.1;
tau = 1*1e-6;
niter = 70000;



p = 2;
tau = 1e-6;
niter = 5000;



p = 3;
tau = 1e-8;
niter = 100000;

p = 1;
tau = 1*1e-6;
niter = 70000;

q = 100;
it_disp = round(linspace(1,niter,q));
idisp = 1;
f = f0; F = f;
for it=1:niter
    if it==it_disp(idisp)        
        clf; hold on;
        for jt=1:idisp
            s = (jt-1)/(q-1);
            col = [s 0 1-s];
            plot(x,F(:,jt), 'color', .5 + .5*col, 'LineWidth', 1);
        end    
        plot(x,f, 'color', col, 'LineWidth', 3);
        axis([0 1 0 1.01]);
        box on;
        set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
        drawnow;
        mysaveas(idisp); 
        idisp = idisp+1;
        F(:,end+1) = f;
    end
    %
    f = f + tau*lapl(f,p);
end
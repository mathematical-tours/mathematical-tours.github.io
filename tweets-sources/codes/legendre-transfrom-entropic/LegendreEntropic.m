%% 
% Legendre transform via hopf-cole log-transform

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 1024; 

x = linspace(-1,1,n)';
epsilon = (.2)^2;
K0 = @(epsilon)exp(-(x-x').^2/(2*epsilon));
K = @(epsilon)diag(1./sum(K0(epsilon),2)) * K0(epsilon);
C = @(f,epsilon)exp( (x.^2/2-f)/epsilon );
Ci = @(F,epsilon)x.^2/2 - epsilon*log(F);

u = 1;
f = x.^2 + .3*u*sin(6*pi*x);

epsilon = (.021).^2;
clf; hold on;
plot(x, f, 'k');
plot(x, Lgd(Lgd(f,epsilon),epsilon), 'b');

Lgd = @(f,epsilon)Ci( 1 ./ (K(epsilon)*C(f,epsilon)), epsilon );

q = 40;
eps_list = ( 3*linspace(1,0,q).^3 + .021 ) .^2;
clf; hold on;
plot(x, f, 'k', 'LineWidth',2);
for it=1:q
    s = (it-1)/(q-1);
    epsilon = eps_list(it);
    plot(x, Lgd(Lgd(f,epsilon),epsilon), 'color', [s 0 1-s]);
    axis([-1, 1, min(f)-.1, max(f)+.1]);
    set(gca, 'XTick', [], 'YTick', []);
    box on;
    drawnow;
    mysaveas(it);
end

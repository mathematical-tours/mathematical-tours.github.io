%%
% test for marcenko-pastur law

rep = ['results/marcenko-pastur/'];

Q = 1000;
T = linspace(1e-5,3,Q);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

K = 10;
blist = linspace(0,1,K);
blist = 1./(1:K);
blist = 1./[1 2 4 10 50 100];
K = length(blist);

lgd = {};
clf; hold on;
for k=1:K
    beta = blist(k);
    col = (k-1)/(K-1);
    a = (1-sqrt(beta))^2;
    b = (1+sqrt(beta))^2;
    MP = sqrt((b-T).*(T-a))./( 2*pi*beta*T );
    MP( T<a | T>b ) = 0;
    MP = Q*MP/sum(MP);
    plot(T,MP, 'LineWidth', 2, 'Color', [col 0 1-col]);
    lgd{end+1} = ['\beta=1/' num2str(1/beta)];
end
axis tight;
legend(lgd);
SetAR(1/2);
axis([0 max(T) 0 10]);
box on;
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'mp-laws.eps'], 'epsc');
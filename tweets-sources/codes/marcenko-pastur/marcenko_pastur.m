%%
% test for marcenko-pastur law


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = ['results/marcenko-pastur/'];
[~,~] = mkdir(rep);


r0 = 100000*10; % #repplications

if not(exist('beta'))
beta = 1/2;
beta = 1/4;
beta = 1/10;
beta = 1/50;
end

slist = [2 5 10 100];


L1 = {};
for is=1:length(slist)
    s = slist(is);
    p = s/beta;
    
    r = round(r0/s);
    S = @(X)X'*X;
    L = [];
    for i=1:r
        progressbar(i,r);
        L = [L; eig(S(randn(p,s))/p)];
    end
    L1{is} = L;
end

a = (1-sqrt(beta))^2;
b = (1+sqrt(beta))^2;

q = 80; Q = 1000;
t = linspace(1e-5,b*1.5,q);
T = linspace(1e-5,b*1.5,Q);


% true limit law
MP = sqrt((b-T).*(T-a))./( 2*pi*beta*T );
MP( T<a | T>b ) = 0;
MP = Q*MP/sum(MP);

lgd = {};
clf; hold on;
for is=1:length(slist)
    col = (is-1)/(length(slist)-1);
    L = L1{is};
    %
    h = hist(L, t);
    plot(t, q*h/sum(h), 'LineWidth', 2, 'Color', [col 0 1-col]);
    lgd{end+1} = ['s=' num2str(slist(is))];
end
plot(T,MP, 'k', 'LineWidth', 3);
box on;
legend(lgd);
set(gca, 'FontSize', 20);
axis tight;
SetAR(1/2);
saveas(gcf, [rep 'marcenko-pastur-' num2str(1/beta) '.eps'], 'epsc');
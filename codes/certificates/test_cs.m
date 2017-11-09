%%
% test for certificate when using Gaussian matrices.


rep = 'results/certificates-cs/';
[~,~] = mkdir(rep);

n = 64;

I = round([.2 .5 .8]*n);
sI = [1 -1 1]';


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

Plist = 1./[1 2 4 8 16]*n;

ms = 20; lw = 2; fs = 15;
for i=1:length(Plist)
    p = Plist(i);
    A = randn(p,n);
    pF = pinv(A(:,I))'*sI;
    etaF = A'*pF;
    clf; hold on;
    plot(etaF, 'k.-', 'MarkerSize', ms, 'LineWidth', lw);
    plot(I, sI, 'r.', 'MarkerSize', 2*ms);
    plot([1 n], [1 1], '--', 'LineWidth', 2);
    plot([1 n],-[1 1], '--', 'LineWidth', 2);
    set(gca, 'FontSize', fs);
    axis([1 n -1.2 1.2]);
    SetAR(1/2); box on;
    saveas(gcf, [rep 'certif-cs-' num2str(p) '.eps'], 'epsc');
end
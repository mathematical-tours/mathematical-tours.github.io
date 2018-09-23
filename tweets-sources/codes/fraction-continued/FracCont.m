%%
% Compute continued fraction expansion

addpath('../toolbox/');
rep = MkResRep();

X = {3/17, sqrt(2), 2^(1/6), pi};
X = {EvalFracCont([0 1 2 3 4 5 6 7 8 9 10]), sqrt(2), 2^(1/6), pi};
X = {(1+sqrt(5))/2, 2^(1/3), pi}; 
X = {pi};
col = distinguishable_colors(length(X));

n = 40;
nmax = 1e7;

KMax = -Inf;
clf; hold on;
for it=1:length(X)
    x = X{it};
    a = ComputeContFrac(x,n);
    
    
    Z = [];
    for i=1:n
        Z(i) = EvalFracCont(a(1:i));
    end
    [H,K] = Convergents(a);
    Z = H./K;
    
    I = find(K<1e15 & log10(abs(Z-x))>-18 );
    KMax = max(KMax,max(K(I)));
    plot(log10(K(I)), log10(abs(Z(I)-x)), '.-', 'Color', col(it,:), 'LineWidth', 2, 'MarkerSize', 25);
    if it==length(X)
        % best approx
        Wi = 1:max(K(I));
        W = round(x*Wi)./Wi;
        plot(log10(Wi), log10(abs(W-x)), '-', 'Color', .5+.5*col(it,:), 'LineWidth', 1);
    end
end
axis([0 log10(KMax) -14 0]);
axis tight;
set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20);
box on;
saveas(gcf, [rep 'FracCont.eps'], 'epsc');

% AutoCrop(rep, 'gauss');
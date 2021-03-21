%%
% Display of Hermite orthogonal basis of L^2


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


% [ 1, 2*x, 4*x^2 - 2, 8*x^3 - 12*x, 16*x^4 - 48*x^2 + 12]

n = 1024;
x = linspace(-10,10,n)';
K = 20; % #hermite functions

syms y;
H = [];
for i=0:K-1
    P = hermiteH(i,y);
    H(:,i+1) = subs(P, x); % plot(x, hermiteH(1,x));
end
Hf = H .* exp(-x.^2/2);
Hf = Hf .* sign(Hf(end/2,:));
Hf = Hf ./ max(Hf);

it = 0;
for i=1:K
    clf; hold on;
    for j=1:i-1
        s = (j-1)/(K-1);
        plot(x,Hf(:,j), 'LineWidth', 2, 'color', .7+.3*[s 0 1-s]);
    end
    s = (i-1)/(K-1);
    plot(x,Hf(:,i), 'LineWidth', 2, 'color', [s 0 1-s]);
    axis([-8 8 -1.3 1.3]);
    box on;
    set(gca,'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    drawnow;
    for m=1:5
        it = it+1;
        mysaveas(it);
    end
end
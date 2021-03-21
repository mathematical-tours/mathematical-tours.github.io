%%
% Display of pseudo spectrum


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 10;
A{1} = randn(n) + 1i*randn(n);
A{2} = randn(n) + 1i*randn(n);

n = 15;
A{1} = gallery('grcar',n);
A{2} = randn(n) + 1i*randn(n);

for k=1:2
    A{k} = A{k}/max(abs(eig(A{k})));
end

m = 150;
t = linspace(-1.05,1.05,m);
Z = t + 1i*t';

q = 80;
for it=1:q
    s = (it-1)/(q-1);
    B = (1-s)*A{1}+s*A{2};
    %
    H = zeros(m);
    for i=1:m*m
        H(i) = norm( inv(B-Z(i)*eye(n)) );
    end
    % display    
    r = 15;
    clf; hold on;
    imagesc(t,t,log10(H));
    contour(t,t,log10(H),linspace(0,2,r), 'k');
    colormap(parula(r-1));
    caxis([0 2]);
    plot(eig(B), 'r.','MarkerSize', 20);
    axis equal; axis([-1 1 -1 1]);
    axis off;
    drawnow;
    mysaveas(it);
end
%%
% Simple legendre transforms.


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep  'a-' name '-' znum2str(it,3) '.png']);


q = 70;

s = linspace(-3,3, 1000);
p = linspace(1+1e-3, 10, q);
p = 1 + 1e-3 + 10 * linspace(0,1, q).^2;
[S,P] = meshgrid(s,p);

Q = P./(P-1); % 1/q = 1-1/p ==> q = P./(P-1);

F = 1./P .* abs(S).^( P );
Fs = 1./Q .* abs(S).^( Q );


% 



for it=1:q
    clf; hold on;
    for jt=1:it-1
        t = (jt-1)/(q-1);
        plot( s,F(jt,:), 'LineWidth', 1, 'Color', .5 + .5*[t 0 1-t] );
    end
    plot( s,F(it,:), 'LineWidth', 3, 'Color', [t 0 1-t] );
    axis([min(s) max(s) 0 3]); box on;
    set(gca, 'XTick', [-1 1], 'YTick', [0 1 2]);
    drawnow;
    %mysaveas('primal', it);
end

for it=1:q
    clf; hold on;
    for jt=1:it-1
        t = (jt-1)/(q-1);
        plot( s,Fs(jt,:), 'LineWidth', 1, 'Color', .5 + .5*[t 0 1-t] );
    end
    plot( s,Fs(it,:), 'LineWidth', 3, 'Color', [t 0 1-t] );
    axis([min(s) max(s) 0 3]); box on;
    set(gca, 'XTick', [-1 1], 'YTick', [0 1 2]);
    drawnow;
   % mysaveas('dual', it);
end


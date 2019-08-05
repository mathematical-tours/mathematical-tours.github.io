%%
% Plot influence of prior on posterior


addpath('../toolbox/');
rep = MkResRep();


A = 4; 
n = 1001;
t = linspace(-2,6,n)';

G = @(m,s)1/sqrt(2*pi*s^2) * exp(-(t-m).^2/(2*s^2));

%    plot( t, LikeL, 'r', 'Linewidth', 2 );
myplot = @(y,c)area( t, y, 'FaceColor', c, 'EdgeColor', c, 'LineWidth', 2);


q = 50; 
for it=1:q
    s = (it-1)/(q-1);
    clf; hold on;
    % likelihood
    LikeL = G(0,.8); 
%    plot( t, LikeL, 'r', 'Linewidth', 2 );
    h = myplot(LikeL,  'r');  h.FaceAlpha = 0.5;
    [vj,j] = max(LikeL);
    plot(t(j),vj, 'r.', 'MarkerSize',40);
    plot([t(j) t(j)],[0 vj], 'r:', 'LineWidth',4);    
    % Prior
    
    m = 3 + cos(s*2*pi);
    sigma = .5 + .2*cos(s*2*pi);
    
    
    m = 3 + cos(s*2*pi);
    sigma = .5 + .2*cos(s*2*pi+pi);
    
    
    Prior = G(m,sigma); 
    h = myplot(Prior, 'b' ); h.FaceAlpha = 0.5;
    % Posterior
    Post = LikeL .* Prior;
    Post = Post / (mean(Post)*(max(t)-min(t)));
    set(gcf, 'alpha', .5);
    h = myplot( Post,  [0 .7 0] ); h.FaceAlpha = 0.5;
    [vj,j] = max(Post);
    plot(t(j),vj, '.', 'MarkerSize',40, 'color', [0 .7 0]);
    plot([t(j) t(j)],[0 vj], ':', 'LineWidth',4, 'color', [0 .7 0]);
    %
    axis([min(t) max(t) 0 1.5]);
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png'] );
end
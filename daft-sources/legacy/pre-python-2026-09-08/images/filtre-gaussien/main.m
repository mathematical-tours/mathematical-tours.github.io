[im, cm] = imread('maurice.png');
colormap(cm);

[n,m] = size(im);

co = [0.001,0.005,0.01,0.02];
nb = length(co);

colormap gray;
    
for s=1:nb
    % le filtre gaussien
    subplot(1,4,s);
    f = compute_filter(n,co(s));
    m = max(max(f));
    p = 50;
    image( f(p:(n-p),p:(n-p))*255/m );
    title( sprintf('t=%.3f', co(s)) );
    if s==1
        ylabel('Fonction g'); 
    end
    colormap gray;
    axis image;
    set(gca, 'XTick', []); set(gca, 'YTick', []);
end

saveas(gcf, '../../images/filtre-gaussien-1', 'eps');
saveas(gcf, '../../images/filtre-gaussien-1', 'png');

pause;

for s=1:nb
    % l'image
    subplot(1,4,s);
    f = compute_filter(n,co(s));
    y = filter2(f,im);
    image(y);
    if s==1
        ylabel('Fonction f*g'); 
    end
    colormap(cm);
    axis image; 
    set(gca, 'XTick', []); set(gca, 'YTick', []);
end 

saveas(gcf, '../../images/filtre-gaussien-2', 'eps');
saveas(gcf, '../../images/filtre-gaussien-2', 'png');

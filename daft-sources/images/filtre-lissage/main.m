clf;
n = 128;


co = [0.00001,0.005,0.01,0.02];
nb = length(co);


h = n/4;
x = ( 0:1/(h-1):1 )';
y = [ones(h,1)+rand(h,1)*0.4; sin(x*pi); ones(h,1)+rand(h,1); x.^2];

% tout d'abord un dirac.
for s=1:nb
    % le filtre gaussien
    subplot(1,4,s);
    f = compute_filter(n+1,co(s));
    p = 40;
    sel = p:(n-p+1);
    x = ( -(n/2-p) ):( n/2-p+1 );
    plot( x, f(sel), 'k.:' );
    axis tight;
    axis square;
    title( sprintf('t=%.3f', co(s)) );
    if s==1
        ylabel('Fonction g'); 
    end
end

saveas(gcf, '../../images/filtre-lissage-1', 'eps');
saveas(gcf, '../../images/filtre-lissage-1', 'png');

pause;

x = 0:n-1;
for s=1:nb
    % l'image
    subplot(1,4,s);
    f = compute_filter(n,co(s));
    f = f( [(n/2+1):n,1:n/2] );
    yy = filter(f,1,y);
    plot( x, yy, 'k.:' );
    axis( [0,n-1,0,1] );
    axis square;
    if s==1
        ylabel('Fonction f*g'); 
    end
end 

saveas(gcf, '../../images/filtre-lissage-2', 'eps');
saveas(gcf, '../../images/filtre-lissage-2', 'png');
% Transformée de Hartley généralisée

n = 16;
p = 512;    % taille avec zero pad.
f = [(0:n/2-1)';(n/2:-1:1)'];

f = [f;zeros(p-n,1)];

nb = 9;
x = (-p/2+1:p/2)*(n/p);
for i=0:nb-1
    lambda = i*2*pi/nb;
    y = htgen(f,lambda);
    y = [y(p/2+1:p);y(1:p/2)];
    subplot(3,3,i+1);
    plot(x, y, 'k'); 
	text( 2,40, ['\', sprintf( 'lambda=%0.1f', lambda )] );   
    axis square;
%    axis tight;
    axis ([-8,8,-60,60]);
end

saveas(gcf, '../transformee-hartley-generalisee', 'eps');
saveas(gcf, '../transformee-hartley-generalisee', 'png');
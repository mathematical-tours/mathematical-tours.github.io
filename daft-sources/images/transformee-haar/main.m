n = 32;

xx = 0:1/(n-1):1;

x = (0:n-1)';
y = haar(x);
subplot(2,2,1);
plot(xx,x,'k:.');
title('Signal original');
subplot(2,2,2);
plot(xx,y,'k:.');
title('Transformée de Haar discrète');

x = sin(xx*2*pi)';
y = haar(x);
subplot(2,2,3);
plot(xx,x,'k:.');
title('Signal original');
subplot(2,2,4);
plot(xx,y,'k:.');
title('Transformée de Haar discrète');

saveas(gcf, '../transformee-haar', 'eps');
saveas(gcf, '../transformee-haar', 'png');
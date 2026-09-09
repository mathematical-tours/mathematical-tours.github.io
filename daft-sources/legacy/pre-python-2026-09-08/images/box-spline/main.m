P = 1024;
x = 4*(-P/2:P/2)/P;
f1 = (x<=0.5)&(x>=-0.5);

sel = [(P/2+1):P+1, 1:P/2];
f = f1(sel);

f2 = real( ifft( fft(f).*fft(f) ) )*4/P;
f2 = f2(sel);

f3 = real( ifft( fft(f).*fft(f).*fft(f) ) )*4/P*4/P;
f3 = f3(sel);

f4 = real( ifft( fft(f).*fft(f).*fft(f).*fft(f) ) )*4/P*4/P*4/P;
f4 = f4(sel);

subplot(2,2,1);
plot(x,f1,'k');
axis([-2,2,0,1.05]);
title('\beta^0');

subplot(2,2,2);
plot(x,f2,'k');
axis([-2,2,0,1.05]);
title('\beta^1');

subplot(2,2,3);
plot(x,f3,'k');
axis([-2,2,0,1.05]);
title('\beta^2');

subplot(2,2,4);
plot(x,f4,'k');
axis([-2,2,0,1.05]);
title('\beta^3');

saveas(gcf, '../../images/box-spline', 'eps')
saveas(gcf, '../../images/box-spline', 'png')
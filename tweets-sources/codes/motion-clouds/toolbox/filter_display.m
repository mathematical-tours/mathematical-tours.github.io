function filter_display(H)

g = ifft(H,[],3);
e = squeeze( sum(sum(abs(g).^2, 1), 2) );
E = sum(H, 3);
clf; 
subplot(2,2,3); plot(fftshift(e)); axis tight;
title('t, time');
subplot(2,2,1); imageplot(fftshift(E));
h = real( ifftn(H) );
title('2D, Fourier');
subplot(2,2,2); imageplot(fftshift(sum(abs(h),3)));
title('2D, Space');
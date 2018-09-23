%%
% test for sinc using fft

n = 64; 
q = 1*30; % oversampling ratio

x = [0:n*q/2-1, -n*q/2+1:-1]/n;

f = real( ifft(abs(x)<=1/3) );

clf;
plot(fftshift(f), '.-', 'LineWidth', 2, 'MarkerSize', 20);
 axis tight;
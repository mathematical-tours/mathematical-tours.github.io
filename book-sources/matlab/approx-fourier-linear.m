F = fftshift(fft2(f)); 
n = size(f,1);
F1 = zeros(n);
sel = (n/2-m:n/2+m)+1;
F1(sel,sel) = F(sel,sel); 
f1 = real( ifft2(fftshift(F1)) );
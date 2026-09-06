F = fft2(f);
a = sort(abs(F(:))); a = a(end:-1:1);
T = a(M+1);
F = F .* (abs(F)>T);
f1 = real( ifft2(F) );
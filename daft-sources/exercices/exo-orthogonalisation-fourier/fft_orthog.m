function y = fft_orthog(x)
y = fft(x);
for i=1:length(y)
    if y(i)==0
        error('La TFD de x s''annule.');
        return;
    end
end
y = y./abs(y);
y = real( ifft(y) );
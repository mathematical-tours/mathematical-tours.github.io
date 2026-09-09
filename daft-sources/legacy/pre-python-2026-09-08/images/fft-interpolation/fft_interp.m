function [x,y] = fft_interp(f,facteur)

[N,NN] = size(f);
NN = N*facteur;
h = floor(N/2);

f_fft  = fft(f);
ff_fft = zeros(NN,1);
ff_fft(1:h+1) = f_fft(1:h+1);
ff_fft((NN-h+1):NN) = f_fft((N-h+1):N);

y = facteur*real(ifft(ff_fft));
x = (0:(NN-1))/facteur;

% on extrait seulement la partie qui nous intéresse
extr = (1:(facteur*(N-1)+1));
y = y(extr);
x = x(extr);
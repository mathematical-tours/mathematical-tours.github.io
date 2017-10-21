%%
% Turn any vector into eigenvector of fft

N = 2048;
x = [0:N/2,-N/2+1:-1]';


%
t = linspace(-1,1,N);
t = x/N;
f = double( abs(t)<.2 );



k = 1;
f = exp( 2i*pi/N*k*(0:N-1)' );


%
s = .1;
f = exp(-t.^2/(2*s^2));

mfft = @(x)fft(x)/sqrt(N);
eigenfix = @(f)( f + mfft(f) + mfft(mfft(f))+ mfft(mfft(mfft(f))) )/4;

plot(eigenfix(f), '.-');

clf;
subplot(2,1,1);
plot(real(eigenfix(f)) );
subplot(2,1,2);
plot(imag(eigenfix(f)) );

norm( mfft(eigenfix(f)) - eigenfix(f) ) / norm( eigenfix(f) )
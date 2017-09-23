%%
% test speed of fft

N = 1024*32;

A = exp( -2i*pi/N * (0:N-1)'*(0:N-1) );

clf;
imagesc(real(A));
imagesc(real(fftshift(A)));

u = randn(N,1);
tic;
v = A*u;
a = toc;
%
tic;
w = fft(u);
b = toc;

norm(v-w)/norm(v)
fprintf('Time gain: %.3f\n', a/b );
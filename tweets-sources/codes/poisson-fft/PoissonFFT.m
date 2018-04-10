%%
% Fast Poisson solver.

n = 128;

x = randn(2*n,1);
y = fft(x);
y = imag(y)*1i; y(1:2:end) = 0;
x1 = ifft(y);

% RHS
a = .3; b = -.4;
y = zeros(2*n,1);

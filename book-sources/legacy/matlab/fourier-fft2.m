% Compute the Fourier transform.
F = fft2(f);  n = size(f,1);
% Compute marsked Fourier transform.
t = linspace(-pi(),pi(),n);
h = (cos(t)+1)/2; h = h'*h;
F1 = fft2(f.*h);
% Compute Log of Fourier transforms.
L = fftshift(log( abs(F)+1e-1 ));
L1 = fftshift(log( abs(F1)+1e-1 ));
% display
clf; imageplot( {L L1}, {'FFT' 'Masked FFT'} );

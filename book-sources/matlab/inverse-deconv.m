% Fourier transform of the filter, assumed to be symetric
hF = fft2(h);
% Shortcut for the filtering operator
filt = @(x)real(ifft2(fft2(x).*hF));
% Iterative soft thresholding in wavelets.
fspars = y; % initialization
tau = 1.5/max(abs(hF(:)));
for i=1:niter
    fspars = fspars + tau * filt( y-filt(fspars) );
    fW = perform_wavelet_transf(fspars, Jmin, +1,options);
    fW = perform_thresholding( fW, lambda*tau, 'soft' );
    fspars = perform_wavelet_transf(fW, Jmin, -1,options);
end
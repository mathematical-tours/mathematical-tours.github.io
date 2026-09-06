% Fourier transform of the filter
phiF = fft2(phi);
% Perform the L2 inversion.
fl2 = real( ifft2( yF .* conj(phiF) ./ ( abs(phiF).^2 + epsilon) ) );
% Compute the Sobolev prior penalty (rescale to [0,1]).
x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
\rho = (X.^2 + Y.^2)*(2/n)^2;
% Perform the sobolev inversion.
fsob = real( ifft2( fft2(y) .* conj(phiF) ./ ( abs(phiF).^2 + lambda*\rho) ) );
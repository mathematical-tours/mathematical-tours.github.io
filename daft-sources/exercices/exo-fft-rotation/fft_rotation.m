function y = fft_rotation(x, theta)
n = size(x);
y = fft_transvec_x( x,-tan(theta/2) );
y = fft_transvec_y( y,sin(theta) );
y = fft_transvec_x( y,-tan(theta/2) );
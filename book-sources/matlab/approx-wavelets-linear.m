fw = perform_wavelet_transf(f,j0,+1);
n = size(f,1); fw1 = zeros(n);
fw1(1:n/4,1:n/4) = fw(1:m,1:m);
f1 = perform_wavelet_transf(fw1,j0,-1);
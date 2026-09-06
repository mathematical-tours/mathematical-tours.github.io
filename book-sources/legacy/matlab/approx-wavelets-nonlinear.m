fw = perform_wavelet_transf(f,j0,+1);
a = sort(abs(fw(:))); a = a(end:-1:1);
T = a(M+1);
fw1 = fw .* (abs(fw)>T);
f1 = perform_wavelet_transf(fw1,j0,-1);
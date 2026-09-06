% Forward wavelet transform
fw = perform_wavelet_transf(f,j0,+1);
% Hard thresholding.
fw = fw .* (abs(fw)>T);
% Backward wavelet transform
f1 = perform_wavelet_transf(fw,j0,-1);
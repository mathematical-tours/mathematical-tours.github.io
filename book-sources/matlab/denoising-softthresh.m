% Forward wavelet transform
fw = perform_wavelet_transf(f,j0,+1);
% Soft thresholding.
fw1 = max(1-T./(abs(fw)+1e-20), 0).*fw;
% Add back the coarse scale coefficients
fw1(1:2^j0) = fw(1:2^j0);
% Backward wavelet transform
f1 = perform_wavelet_transf(fw,j0,-1);
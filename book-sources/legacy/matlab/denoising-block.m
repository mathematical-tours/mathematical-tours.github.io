n = size(f,1);
% Forward wavelet transform.
fw = perform_wavelet_transf(f,j0,+1);
% Compute indexing of the blocks.
[dX,dY,X,Y] = ndgrid(0:s-1,0:s-1,1:s:n-s+1,1:s:n-s+1);
I = X+dX + (Y+dY-1)*n;
% Extract the blocks
H = reshape(fw(I(:)),size(I));
% Compute the average energy of each block, and duplicate.
v = mean(mean(abs(H).^2,1),2);
v = repmat( max3(v,1e-15), [w w]);
% Stein threshold the blocks.
HT = max3(1-T^2*v.^(-1),0) .* H;
% Reconstruct the thresholded coefficients.
fw(I(:)) = HT(:);
% Inverse wavelet transform.
f1 = perform_wavelet_transf(fw,j0,-1);
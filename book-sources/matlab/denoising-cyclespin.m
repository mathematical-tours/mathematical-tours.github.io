[dY,dX] = meshgrid(0:s-1,0:s-1);
f1 = f*0;
for i=1:m^2
    fs = circshift(f,[dX(i) dY(i)]);
    fw = perform_wavelet_transf(fs,j0,1);
    fw = fw .* (abs(fw)>T);    
    fs = perform_wavelet_transf(fw,j0,-1);
    fs = circshift(fs,-[dX(i) dY(i)]);
    f1 = (i-1)/i*f1 + 1/i*fs;
end
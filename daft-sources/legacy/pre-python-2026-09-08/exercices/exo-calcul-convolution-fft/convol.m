function y = convol(f,g)
n = length(f); m = length(g); p = n/m;
y = zeros(m+n-1,1); 
fg = fft( [g; zeros(m-1,1)] );
for j=0:p-1
   fj = [f( (1:m)+j*m ); zeros(m-1,1)];
   sel = (j*m+1):((j+2)*m-1); % les indices concernés
   y(sel) = y(sel) + ifft( fft(fj).*fg );
end
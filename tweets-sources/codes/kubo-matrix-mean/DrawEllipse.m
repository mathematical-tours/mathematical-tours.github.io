function DrawEllipse(A,m,col,al)

[U,D] = svd(A); D = diag(D);
u = U(1,1) + U(2,1)*1i;
v = U(1,2) + U(2,2)*1i;
%
t = linspace(0,2*pi,128);
E = m + D(1)*cos(t)*u + D(2)*sin(t)*v;
fill(real(E), imag(E), col, 'FaceAlpha',al, 'LineStyle','None');

end
function h = parzen2d(x,y,p,s)

% (x,y): positions in [0,1)
% p: #pixels
% s: width in pixels

xI = round(x*p)+1; yI = round(y*p)+1;
I = find(xI>0 & xI<=p & yI>0 & yI<=p);
xI = xI(I); yI = yI(I);
h = full( sparse(xI,yI,ones(length(xI),1),p,p) );

cconv = @(x,y)real(ifft2(fft2(x).*fft2(y)));

if 1
    u = [0:p/2, -p/2+1:-1]';
    g = exp(-u.^2/(2*s^2));
    g = g*g'; h = cconv(h,g);
else
    u = linspace(-p/2,p/2,p);
    g = exp(-u.^2/(2*s^2));
    h = conv2(conv2(h,g,'same'),g','same');
    h = h/sum(h(:));
end

end
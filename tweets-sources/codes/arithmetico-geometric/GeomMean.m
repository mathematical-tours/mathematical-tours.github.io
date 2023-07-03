function m = GeomMean(x,y)

a = angle(x);
b = angle(y);
if a-b>pi
    a = a-2*pi;
end
if a-b<-pi
    a = a+2*pi;
end
if abs(a-b)>pi
    warning(problem);
end
m = sqrt(abs(x*y)) * exp( 1i*(a+b)/2 );
end
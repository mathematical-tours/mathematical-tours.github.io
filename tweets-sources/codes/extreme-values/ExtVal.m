function h = ExtVal(s,xi)

n = length(s);
A = max(s);

if not(xi==0)
    f = exp(-(1+xi*s).^(-1/xi)) .* ((xi*s+1)>0);
else
    f = exp(-exp(-s));
end
h = diff(f); h = (n-1)/(2*A)*h; % /sum(h);
h = max(real(h),0);
end
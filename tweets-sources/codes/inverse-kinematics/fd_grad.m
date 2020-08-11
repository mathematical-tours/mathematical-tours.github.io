function dt = fd_grad(f,t)

n = length(t);
delta = 1e-12;
for i=1:n
    t1 = t; t1(i) = t(i)+delta;
    dt(i) = (f(t1)-f(t))/delta;
end
dt = dt(:);

end
function g = LinApprox(f,p)

n = size(f,1);
g = zeros(n);
q = ceil(n/p);
for i=1:p
    for j=1:p
        si = min((i-1)*q+1:i*q,n);
        sj = min((j-1)*q+1:j*q,n);
        g(si,sj)=mean(mean(f(si,sj)));
    end
end

end
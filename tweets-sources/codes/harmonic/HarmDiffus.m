function f = HarmDiffus(f, x, a, niter, ndisp, ord)

if nargin<2
    ord = 1;
end

n = size(f,1);
for it=1:niter
    
    
    L = f;
    for k=1:ord
        L = L - ( L([2:end end],:) + L([1 1:end-1],:) + L(:,[2:end end]) + L(:,[1 1:end-1]) )/4;
    end
    tau = .5;
    f = f - tau * L;
    
    %
    f(real(x) + (imag(x)-1)*n) = a;
    if mod(it,ndisp)==1
        clf;
        imageplot(f);
        drawnow;
    end
end

end
function f = HarmDiffus(f, x, a, niter, ndisp)

n = size(f,1);
for it=1:niter
    f = ( f([2:end end],:) + f([1 1:end-1],:) + f(:,[2:end end]) + f(:,[1 1:end-1]) )/4;
    f(real(x) + (imag(x)-1)*n) = a;
    if mod(it,ndisp)==1
        clf;
        imageplot(f);
        drawnow;
    end
end

end
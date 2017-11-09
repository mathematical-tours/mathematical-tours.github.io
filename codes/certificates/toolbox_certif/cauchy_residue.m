%computes the contour integral of F, given the poles p inside the contour
function res = cauchy_residue(F,p)
syms z;

res = 0;

for k=1:length(p)
    f = (z-p(k))*F;
    i=0;
    while true
        try
            
             subs(f, z, p(k));   
             break;
        catch
            f = (z-p(k))*f;
             i=i+1;
        end
        
    end
    dfi = diff(f,z,i)/factorial(i);
    res = res + 2*pi*1i*subs(dfi, z, p(k));
end


end
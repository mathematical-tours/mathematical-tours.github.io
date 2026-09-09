function ff = antisymetrise(f)
n = length(f)+1; ff = zeros(2*n,2*n,2*n);
for x=1:2*n
for y=1:2*n
for z=1:2*n
    if mod(x-1,n)==0 | mod(y-1,n)==0 | mod(z-1,n)==0
        ff(x,y,z) = 0;
    else
        signe = 1; nx = x; ny = y; nz = z;   
        if(x>n) signe = -signe; nx = 2*n-x+2; end
        if(y>n) signe = -signe; ny = 2*n-y+2; end
        if(z>n) signe = -signe; nz = 2*n-z+2; end
        ff(x,y,z) = signe*f(nx-1,ny-1,nz-1);
    end
end
end
end
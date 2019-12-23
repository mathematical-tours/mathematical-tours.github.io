function g = ErrorDiffusion(f)

n = size(f,1);

method = 'floyd';

g = f;
for i=2:n-1
    for j=2:n-1
        e = g(i,j) - (g(i,j)>.5);
        g(i,j) = g(i,j)>.5;
        switch method
            case 'simple'
                g(i+1,j) = g(i+1,j) + e/2;
                g(i,j+1) = g(i+1,j) + e/2;
            otherwise
                e = e/16;
                g(i,j+1) = g(i,j+1)     + e*7;
                g(i+1,j+1) = g(i+1,j+1) + e;
                g(i+1,j)   = g(i+1,j)   + e*5;
                g(i+1,j-1) = g(i+1,j-1) + e*3;
        end
    end    
end

end

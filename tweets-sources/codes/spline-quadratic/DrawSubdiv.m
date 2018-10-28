function DrawSubdiv(f,p, Col)

n = length(f);

lw = 1;

hold on;
for i=1:n
    a = f(i);
    b = f(mod(i-2,n)+1); b = (a+b)/2;
    c = f(mod(i,n)+1); c = (a+c)/2;
    for k=1:p
        t = (k-1)/(p-1);
        plot( (1-t)*[a c] + t*[b a], 'Color', Col, 'LineWidth', lw );        
    end
    
end

end
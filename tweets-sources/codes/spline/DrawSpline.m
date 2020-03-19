function DrawSpline(f,p,deg,Col)

f = f(:);
n = length(f);
n = floor(n/(deg-1))*(deg-1);

lw = 2;

hold on;
for i=1:(deg-1):n
    I = mod( (i:i+deg)-1,n )+1;
    fI = f(I);
    fI = [ (fI(1)+fI(2))/2; fI(2:end-1); (fI(end-1)+fI(end))/2 ];
    g = EvalSpline(fI);
    plot(fI, '.-', 'Color', [1 1 1]*.6, 'LineWidth', 2, 'MarkerSize', 25);
    plot(g,'Color', Col, 'LineWidth', lw);
end


end



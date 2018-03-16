function PlotLines(c,cn, str)

q = 100;
R = 20;
t = linspace(-R,R,q);
hold on;
for i=1:length(c)
    plot(real(c(i) + t*cn(i)), imag(c(i) + t*cn(i)), 'color', [1 1 1]*.7);
end

end
function DrawRod(x)

n = length(x);
hold on;
for i=1:n-1
    t = (i-1)/(n-2);
    plot(real(x(i:i+1)),imag(x(i:i+1)), '.-', 'color',[t 0 1-t], 'LineWidth', 2, 'MarkerSize', 25);
end
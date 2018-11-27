function p = DrawPoly()

p = [];
clf; hold on;
while true
    plot(p);
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    p(end+1) = a+1i*b;
end
p = p(:);

end
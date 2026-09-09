N = 256;
x = -1:2/(N-1):1;

for i=1:6
    y = cos( i*acos(x) );
    subplot(2,3,i);
    plot(x,y);
    str = sprintf('T_%d', i);
    title(str);
end

saveas(gcf, '../polynomes-chebyshev', 'eps')
saveas(gcf, '../polynomes-chebyshev', 'png')
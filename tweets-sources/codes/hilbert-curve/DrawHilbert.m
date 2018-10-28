
addpath('../toolbox/');
rep = MkResRep();

r = 1;
for n=1:6
    [x,y] = hilbert(n);
    lw = 4;
    k = length(x);
    clf; hold on;
    for i=1:k-1
        t = (i-1)/(k-2);
        plot(x(i:i+1),y(i:i+1), 'LineWidth', lw, 'Color', [1-t, 0, t]);
    end
    axis equal; axis off;
    for i=1:5
        saveas(gcf, [rep 'hilbert-' znum2str(r,2) '.png']);
        r = r+1;
    end
end
% AutoCrop(rep, 'hilbert-');
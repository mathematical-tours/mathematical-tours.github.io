
if not(exist('test'))
    test = 0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));


clf; hold on;
f0 = [];
while true
    axis equal; axis([0 1 0 1]);  % axis off;
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    f0(end+1) = a + 1i*b;
    plot(f0, 'r', 'LineWidth', 2);
end
k = size(f0,2);
plot(f0([1:end,1]), 'r', 'LineWidth', 2);
f1 = [];
for it=1:k
    axis([0 1 0 1]); axis equal; axis off;
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    f1(end+1) = a + 1i * b;
end
f0 = f0(:); f1 = f1(:);

p = 10;

q = 50;  % animation
for deg=2:3
for it=1:q
    t = (it-1)/(q-1);
    ft = f0*(1-t) + f1*t;
    
    clf;
    DrawSpline(ft,p,deg,[t 0 1-t]);
    axis equal; axis([0 1 0 1]); set(gca, 'Xtick', [], 'Ytick', []);
    drawnow;
    saveas(gcf, [rep  'anim-' num2str(deg) '-' znum2str(it,2) '.png'], 'png');
end
end

% AutoCrop(rep,'anim');


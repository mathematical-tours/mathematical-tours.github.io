function [f0,f1,f2] = selectPoly()


%%
% First we create a 2D closes polygon.

clf; hold on;
f0 = [];
while true
    axis equal; axis([0 1 0 1]);  % axis off;
    box on; set(gca,'XTick',[],'YTick',[]);
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    f0(end+1) = a + 1i*b;
    plot(f0, 'r', 'LineWidth', 2);
end
%%
clf; hold on;
k = size(f0,2);
plot(f0([1:end,1]), 'r', 'LineWidth', 2);
f1 = [];
for it=1:k
    axis equal; axis([0 1 0 1]);  % axis off;
    box on; set(gca,'XTick',[],'YTick',[]);
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    f1(end+1) = a + 1i * b;
end
%%
clf; hold on;
k = size(f0,2);
plot(f0([1:end,1]), 'r', 'LineWidth', 2);
f2 = [];
for it=1:k
    axis equal; axis([0 1 0 1]);  % axis off;
    box on; set(gca,'XTick',[],'YTick',[]);
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    f2(end+1) = a + 1i * b;
end
f0 = f0(:); f1 = f1(:); f2 = f2(:);

end
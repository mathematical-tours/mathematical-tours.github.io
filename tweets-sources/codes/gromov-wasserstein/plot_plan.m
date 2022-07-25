function plot_plan(x,y, gamma)

n = length(x);
m = length(y);
a = ones(n,1)/n;
b = ones(m,1)/m;

t = (0:n-1)/n;
c = sin(pi*t).^2;
colx = [1 0 0]' * c + [0 0 1]' * (1-c);
% barycentric projection from a to b
P = gamma' ./ b(:);
coly = ( P * colx' )';
Q = gamma ./ a(:);
%
hold on;
plotcol(x,colx, 2);
plotcol(y,coly, 2);
u = 60;
I = round(linspace(0,n,u+1)'); I(1)=[];
y1 = Q*y;
for k=1:u
    s = (k-1)/(m-1);
    i = I(k);        
    tau = .1;
    plot(x(i), '.', 'color', colx(:,i), 'MarkerSize', 15);
    plot(y1(i), '.', 'color', colx(:,i), 'MarkerSize', 15);
    h = plot([x(i) y1(i)], '-', 'color', colx(:,i), 'LineWidth', 1);
    h.Color(4) = 0.25;
end
e = .05;
axis equal;  axis([-e,1+e,-e,1+e]);
set(gca, 'XTick', [], 'YTick', []); box on;

end
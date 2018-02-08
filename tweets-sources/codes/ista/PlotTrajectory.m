function PlotTrajectory(x,y, lw)

if nargin<3
    lw = 2;
end
n = length(x);
L = ( x(2:end)-x(1:end-1) ).^2 + ( y(2:end)-y(1:end-1) ).^2;
L = sqrt(L); L = [0; cumsum(L(:))];
L = L/L(end);

for i=1:n-1
    plot(x(i:i+1),y(i:i+1), '-', 'LineWidth', lw, 'color', [L(i) 0 1-L(i)]);
end

end
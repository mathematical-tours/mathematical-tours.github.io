function PlotTrajectory3D(x,y,z, lw, c0, c1)

if nargin<3
    lw = 2;
end
n = length(x);
L = ( x(2:end)-x(1:end-1) ).^2 + ( y(2:end)-y(1:end-1) ).^2;
L = sqrt(L); L = [0; cumsum(L(:))];
L = L/L(end);

for i=1:n-1
	plot3(x(i:i+1),y(i:i+1),z(i:i+1), '-', 'LineWidth', lw, 'color', (1-L(i))*c0 + L(i)*c1);
end

end
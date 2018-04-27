function PlotMesh2d(F,Y)

m = size(F,2);
r = @(u)reshape(u,[2 m]);

lw = 2; ms = 25;
hold on;
for i=1:3
    j = mod(i,3)+1;
    plot(r(Y(1,F([i j],:))), r(Y(2,F([i j],:))), 'k.-', 'LineWidth', lw, 'MarkerSize', ms);
end

end
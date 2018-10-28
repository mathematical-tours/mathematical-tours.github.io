function PlotFlow(X)

co = distinguishable_colors(size(X,1));
set(groot,'defaultAxesColorOrder',co);


d = abs(X(:,1:end-1)-X(:,2:end));
I = find(d>.5);
X(I) = NaN;

hold on;
plot(transpose(X(1:end/2,:)), '-', 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = 1; % keep same color
h = plot(transpose(X(end/2+1:end,:)), '-', 'LineWidth', 2);

end
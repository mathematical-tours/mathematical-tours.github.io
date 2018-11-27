function DisplayDensity(f,s,X)

n = size(f,1);

r = 10;
m = linspace(0,1,r-1)';
CM = (1-m)*[s 0 1-s] + m*[1 1 1];
CM = (m)*[s 0 1-s] + (1-m)*[1 1 1];
%
f = f/max(f(:)); f= rescale(f);
[x,y] = ind2sub([n,n],X);
clf; hold on;
imagesc(f);
colormap(CM);
contour(f,linspace(0,1,r), 'Color', 'k', 'LineWidth', 1); % [s 0 1-s]
plot(y,x, 'k.', 'MarkerSize', 30);
axis image;  axis off; drawnow;
    
end
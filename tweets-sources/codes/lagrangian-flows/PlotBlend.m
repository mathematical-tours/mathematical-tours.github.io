function PlotBlend(X,col)

m = size(X,3);
n = size(X,2);
hold on;
for i=1:m-1
    t = i/(m-1);
    plot(squeeze(X(1,:,i:i+1))', squeeze(X(2,:,i:i+1))', 'color', col*t+(1-t));
end
S = ones(n,1)*20;
scatter(squeeze(X(1,:,end))', squeeze(X(2,:,end))', S, col, 'filled' );

end
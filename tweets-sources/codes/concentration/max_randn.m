%%
% Max of randn variables

n = 5000;

nrep = 100;

clf; hold on;
for i=1:nrep
    Y = cummax(abs(randn(n,1)));
    plot(1:n,Y, 'r', 'LineWidth', 1/2);
    plot(1:n,-Y, 'r', 'LineWidth', 1/2);
end
X = randn(n,1);
%Y = cummax(abs(X));
stem(1:n,X, 'b.');
%plot(1:n,Y, 'r', 'LineWidth', 2);
%plot(1:n,-Y, 'r', 'LineWidth', 2);
plot(1:n,sqrt(2*log(1:n)), 'k--', 'LineWidth', 2);
plot(1:n,-sqrt(2*log(1:n)), 'k--', 'LineWidth', 2);
box on;
saveas(gcf, 'max-gaussian.png', 'png');

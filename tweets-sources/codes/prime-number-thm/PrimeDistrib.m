% https://en.wikipedia.org/wiki/Prime-counting_function

n = 10000;
P = primes(n);

U = zeros(n,1);
for i=1:n
    U(i) = sum(P<=i);
end    

z = (1:n)';
U0 = z./log(z);
clf;
hold on;
plot( z, U, 'b', 'LineWidth', 2 );
plot( z, U0, 'r', 'LineWidth', 2 );

clf; 
semilogx(z, U./U0, 'LineWidth', 2); hold on;
semilogx(z, U./logint(z), 'LineWidth', 2);
semilogx([1 n], [1 1], 'k--', 'LineWidth', 2);
axis([1 n 0 1.3]); box on;
set(gca, 'FontSize', 20);
saveas(gcf, '../results/tmp.png');
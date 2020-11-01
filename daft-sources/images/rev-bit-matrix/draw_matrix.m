function draw_matrix(p)

N = 2^p;
M = zeros(N,1);

for k=1:N
    M(k) = rev_index(p,k-1)+1;
end

plot(1:N, M, 'k.');
axis([1,N,1,N]);
axis square;
str = sprintf('p=%d',p);
title(str);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

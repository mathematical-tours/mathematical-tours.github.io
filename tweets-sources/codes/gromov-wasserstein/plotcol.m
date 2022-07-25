function plotcol(x,col, lw)

% col should be of size 3 x n

if nargin<3
    lw=1;
end

n = length(x);
hold on;
for i=1:n
    j = mod(i,n)+1;
    c = (col(:,i) + col(:,j))/2;
    plot( x([i,j]), 'color', c, 'LineWidth', lw);
end

end
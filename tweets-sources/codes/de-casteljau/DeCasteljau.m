rep = '../results/de-casteljau/';
[~,~] = mkdir(rep);

if not(exist('cnt'))
    cnt = 0;
end
cnt = cnt+ 1;

% click selection
it = 0;
clf; hold on;
Z = [];
while true
    axis([0 1 0 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(end+1) = a+1i*b;
end
Z = Z(:);

k = length(Z);
n = 256;
t = linspace(0,1,n);
B = repmat(Z, [1 n]);
for i=1:k-1
    B = B(1:end-1,:).*repmat(t, [k-i 1]) + B(2:end,:).*repmat(1-t, [k-i 1]);
end

clf; hold on;
plot(Z, 'k-', 'LineWidth', 1);
plot(Z, '.r', 'MarkerSize', 25);
plot(B, 'b', 'LineWidth', 2);
axis([0 1 0 1]); axis equal; axis off;
saveas(gcf, [rep 'bezier-' num2str(cnt) '.eps'], 'epsc');

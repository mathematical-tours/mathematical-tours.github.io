N = 6;
W = ones(1,1);

for i=1:N
	W = [W,W;W,-W];
end

image((W+1)*128);
colormap(gray);
axis image;
grid on;
set(gca,'GridLineStyle', '-');
set(gca, 'XTick', (1:2^N)+0.5);
set(gca, 'YTick', (1:2^N)+0.5);
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');

saveas(gcf, '../../images/matrices-walsh', 'eps');
saveas(gcf, '../../images/matrices-walsh', 'png');
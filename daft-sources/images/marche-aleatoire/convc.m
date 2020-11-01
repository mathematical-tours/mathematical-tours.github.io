% convole n fois de façon cyclique, puis dessine.
function convc(f,v,n)

for i=1:n
    f = real( ifft( fft(f).*fft(v) ) );
end

plot((1:length(f))', f, '*');
axis([1,length(f),0,1]);

set(gca, 'XTick', 1:length(f));
% set(gca, 'YTick', []);
set(gca,'XTickLabel','');
% set(gca,'YTickLabel','');
clf;
n = 16;
p = 16;
i=0;

for x = 1:-0.7/(p-1):0.3

i=i+1;
subplot(4,p/4,i);

fmm = fm^x;
af = abs(real(fmm));
m = max(max( af ));
image( af/m*256 );
colormap(gray);
axis image;
axis off;
t = sprintf( '%0.2f', x );
text( 30,-10, t );

end

clear;
N = 16;
f = (0:N-1)'*(0:N-1);
fm = exp(2*i*pi/N).^f;

P = 4;
f = [ones(1,P),zeros(1,N-2*P+1),ones(1,P-1)]';

clf;
p = 8;
i=0;

for x = 1:-1/(p-1):0
  i=i+1;
  subplot(4,p/4,i);
  ff = (fm^x)*f;
  plot(abs(ff));
  t = sprintf( '%0.2f', x );
  text( 30,-10, t );
end
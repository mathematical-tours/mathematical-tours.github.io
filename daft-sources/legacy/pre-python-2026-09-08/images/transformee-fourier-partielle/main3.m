clear;
N = 17;
x = (0:N-1)'*(0:N-1);
F = 1/sqrt(N)*exp(2*i*pi/N).^x;

V = zeros(N,N);



for i=0:(N-1)/4-1

x = [0,zeros(1,i),1,zeros(1,N-3-2*i),1,zeros(1,i)]';
y = [0,zeros(1,i),1,zeros(1,N-3-2*i),-1,zeros(1,i)]';
vx1 = x+F*x;
vx2 = x-F*x;
vy1 = y+1i*F*y;
vy2 = y-1i*F*y;
vx1 = vx1/norm(vx1);
vx2 = vx2/norm(vx2);
vy1 = vy1/norm(vy1);
vy2 = vy2/norm(vy2);
V(:,1+4*i)=vx1;
V(:,2+4*i)=vx2;
V(:,3+4*i)=vy1;
V(:,4+4*i)=vy2;

end

i = (N-1)/4;
x = [0,zeros(1,i),1,zeros(1,N-3-2*i),1,zeros(1,i)]';
vx1 = x+F*x;
V(:,N) = vx1;

image(255*abs(V^(-1)*F*V))

D = V^(-1)*F*V;
Vv = ctranspose(V);
Id = Vv*V;
m = max(max(abs(Id)));
image(255*abs(Id)/m);
p = 53;


p0 = 25;
f = zeros(p,1);
f(p0) = 1;
psi = zeros(p,1);
p0 = 25;
psi(1)=2;
psi(2)=-1;
psi(p)=-1;
display_transfo(f,psi,3,1);

p0 = 8;
f = [ones(p0,1); zeros(p-2*p0,1); -ones(p0,1)];
psi = zeros(p,1);
psi(1)=1;
psi(2)=1;
psi(p)=-1;
psi(p-1)=-1;
display_transfo(f,psi,3,2);

f = cos( (0:p-1)/(p)*6*pi );
psi = zeros(p,1);
psi(1)=2;
psi(2)=-1.2;
psi(p)=-1.2;
psi(3)=0.2;
psi(p-1)=0.2;
display_transfo(f,psi,3,3);

saveas(gcf, '../ondelettes-corps-finis', 'eps');
saveas(gcf, '../ondelettes-corps-finis', 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saisie des points
axes('position', [0,0,1,1]); 
% [x,y] = ginput;
x= [0.6606,
    0.5345,
    0.3048,
    0.2306,
    0.1995,
    0.2150,
    0.3187,
    0.4482,
    0.7228,
    0.7884,
    0.7798,
    0.6641,
    0.4914,
    0.4378,
    0.4655,
    0.5881,
    0.7280,
    0.7003,
    0.5915];
y =[0.0892,
    0.0817,
    0.1621,
    0.2952,
    0.5214,
    0.7274,
    0.8078,
    0.8656,
    0.8555,
    0.7224,
    0.5992,
    0.6168,
    0.6269,
    0.5013,
    0.3379,
    0.2676,
    0.1420,
    0.0867,
    0.0716];
[N,NN] = size(x);

% un polygone légèrement modifié et tourné
theta = pi*0.2;
noise = 0.1;
xx = x + noise*rand(N,1);
yy = y + noise*rand(N,1);
cpx = complex(xx,yy);
g = mean(cpx);
cpx = g + (cpx-g)*exp(i*theta);
xx = real(cpx);
yy = imag(cpx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dessin des polygones
subplot(1,2,1);
plot([x;x(1)],[y;y(1)],'k', [xx;xx(1)],[yy;yy(1)],'k:');
legend('Forme originale', 'Forme modifiée');
title('Polygones');
axis square
%axis('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des fourier descriptors
fd  = comp_fd(x,y);
fd2 = comp_fd(xx,yy);
[N,NN] = size(fd);
subplot(1,2,2);
plot(1:N,abs(fd),'k.',1:N,abs(fd2),'kx');
axis([1,N,0,2.5]);
axis square
legend('Forme originale', 'Forme modifiée');
title('Descripteurs de Fourier');

plot(real(fd),imag(fd),'k.',real(fd2),imag(fd),'kx');

saveas(gcf, '../../images/fourier-descriptors', 'eps');
saveas(gcf, '../../images/fourier-descriptors', 'png');



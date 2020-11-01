CLF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saisie des points
axes('position', [0,0,1,1]); 
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
[x,y] = ginput;
[N,NN] = size(x);
P = complex(x,y);
nb_iter = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtre qui moyenne
f = zeros(N,1);
f(1) = 0.5;
f(2) = 0.5;
PP = iter_filtre(P,f, nb_iter);
subplot(1,2,1);
plot(PP,'k');
axis square;
axis tight;
title('Filtre moyenne');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtre qui transforme en cercle
f = ones(N,1).*0.8;
f(2) = 1;
f = ifft(f);
PP = iter_filtre(P,f, nb_iter);
subplot(1,2,2);
plot(PP,'k');
axis square;
axis tight;
title('Filtre laissant passer une fréquence');


saveas(gcf, '../../images/filtre-polygones-2', 'eps');
saveas(gcf, '../../images/filtre-polygones-2', 'png');

pause;
CLF;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtre qui explose
% nb_iter = 15;
f = ones(N,1).*0.8;
f(1) = 1.05*exp( complex(0,pi/15) );
f(2) = 1.1*exp( complex(0,pi/10) );
f(3) = 1.05*exp( complex(0,pi/8) );
f(4) = 1.15*exp( complex(0,pi/12) );
f = ifft(f);
PP = iter_filtre(P,f, nb_iter);
subplot(2,2,1);
plot(PP,'k');
axis square;
axis tight;
title('(a) Filtre qui explose');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtre qui explose
% nb_iter = 15;
f = ones(N,1).*0.8;
f(1) = 1.2*exp( complex(0,pi/10) );
f(2) = 1.1*exp( complex(0,pi/8) );
f(3) = 1.15*exp( complex(0,pi/4) );
f(4) = 1.15*exp( complex(0,pi/6) );
f = ifft(f);
PP = iter_filtre(P,f, nb_iter);
subplot(2,2,2);
plot(PP,'k');
axis square;
axis tight;
title('(b) Filtre qui explose');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtre qui converge
nb_iter = 30;
f = ones(N,1).*0.8;
f(1) = 1;
f(2) = 0.95*exp( complex(0,pi/20) );;
f(3) = 0.95;
f = ifft(f);
PP = iter_filtre(P,f, nb_iter);
subplot(2,2,3);
plot(PP,'k');
axis square;
axis tight;
title('(c) Filtre qui converge');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtre qui tourne
% nb_iter = 15;
f = ones(N,1).*0.8;
f(1) = exp( complex(0,pi/15) );
f(2) = exp( complex(0,pi/8) );
f(3) = exp( complex(0,pi/6) );
f = ifft(f);
PP = iter_filtre(P,f, nb_iter);
subplot(2,2,4);
plot(PP,'k');
axis square;
axis tight;
title('(d) Filtre qui tourne');


saveas(gcf, '../../images/filtre-polygones-1', 'eps');
saveas(gcf, '../../images/filtre-polygones-1', 'png');



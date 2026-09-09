% taille du cercle : IMPAIR !
n = 15;


f = zeros(n,1);
f(1) = 1;

nplots = 4;

kval = [0,1,5,20];

% gauche %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = zeros(n,1);
v(2) = 1/2;
v(n) = 1/2;

subplot(nplots, 2, 1);
convc( f,v,kval(1) ); 
text( 12,0.7, sprintf('k=%d',kval(1)) );

subplot(nplots, 2, 3);
convc( f,v,kval(2) ); 
text( 12,0.7, sprintf('k=%d',kval(2)) );


subplot(nplots, 2, 5);
convc( f,v,kval(3) ); 
text( 12,0.7, sprintf('k=%d',kval(3)) );


subplot(nplots, 2, 7);
convc( f,v,kval(4) ); 
text( 12,0.7, sprintf('k=%d',kval(4)) );

% gauche %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = zeros(n,1);
v(1) = 0.3;
v(2) = 0.3;
v(3) = 0.2;
v(n-1) = 0.1;
v(n) = 0.1;

subplot(nplots, 2, 2);
convc( f,v,kval(1) ); 
text( 12,0.7, sprintf('k=%d',kval(1)) );

subplot(nplots, 2, 4);
convc( f,v,kval(2) ); 
text( 12,0.7, sprintf('k=%d',kval(2)) );


subplot(nplots, 2, 6);
convc( f,v,kval(3) ); 
text( 12,0.7, sprintf('k=%d',kval(3)) );


subplot(nplots, 2, 8);
convc( f,v,kval(4) ); 
text( 12,0.7, sprintf('k=%d',kval(4)) );




saveas(gcf, '../marche-aleatoire', 'eps');
saveas(gcf, '../marche-aleatoire', 'png');
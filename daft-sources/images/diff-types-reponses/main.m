% montre les différents types de réponses
clf;
n = 50;

s = 0.05;
x = -1:2/(n-1):1;
f = exp( -(x.^2)/(2*s) );
draw_graph(f, 3, 0, 'exp( - 10 x^2 )');

f = cos( x*pi/2 ).^2;
draw_graph(f, 3, 1, 'cos( 0.5 \pi x )^2' );

f = cos( x*pi/2 );
draw_graph(f, 3, 2, 'cos( 0.5 \pi x )');

% sauvegarde l'image
saveas(gcf, '../diff-types-reponses', 'eps')
saveas(gcf, '../diff-types-reponses', 'png')

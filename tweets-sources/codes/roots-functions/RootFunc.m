% for vizualization
n = 1024;
z = linspace(0,1,n)';

% roots 
p = 200;
x = linspace(.3,.7,p);

F = exp( 1/p * sum(log( abs(x-z)+.001 ),2) );

plot(z, abs(F))
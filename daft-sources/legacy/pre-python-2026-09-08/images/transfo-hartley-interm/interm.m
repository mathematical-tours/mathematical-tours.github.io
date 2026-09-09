function y = interm(l)

global N;
global u1;
global u2;

a = sqrt(N)^l;
y = a*( u1 + (-1)^l*u2 );
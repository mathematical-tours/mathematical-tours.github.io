p = 13;
g = 6;  % c'est bien un générateur de Fp^*, comme le montre le vecteur suivant :
test = mod( g.^(0:p-2) ,p);
x = rand(p,1);
y1 = fft(x);
y2 = fft_chirp(x,g);
norm(y1-y2)
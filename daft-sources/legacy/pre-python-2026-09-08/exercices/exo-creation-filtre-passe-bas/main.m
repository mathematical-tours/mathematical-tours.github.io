% dessine la réponse fréquentielle du filtre (question 2)
N = 64; P = 1024;
f = filtre_passe_bas(N);
ff = [f(1:N/2); zeros(P-N,1); f((N/2+1):N)];
ff = real(fft(ff)); 
plot(ff); axis tight;
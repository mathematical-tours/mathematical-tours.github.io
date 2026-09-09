function f = filtre_passe_bas(N)
f = (0:N-1)'; f = (f<=N/4)|(f>=3*N/4);
f = real(ifft(f));
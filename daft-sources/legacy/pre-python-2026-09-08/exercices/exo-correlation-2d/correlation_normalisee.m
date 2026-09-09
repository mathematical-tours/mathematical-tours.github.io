function y = correlation_normalisee(f,g)
N = length(f); P = length(g);
% renormalisation
g = g - mean(mean(g)); g = g/norm(g);
% calcul du numérateur
ff = zeros(N+P-1,N+P-1); ff(1:N,1:N) = f;
gg = zeros(N+P-1,N+P-1); gg(1:P,1:P) = g;
fg = real(ifft2( fft2(ff).*conj(fft2(gg)) ));
fg = fg(1:N,1:N);
% calcul du dénominateur
s1 = somme_glissante(f,P,1);
s2 = somme_glissante(f,P,2);
denom = sqrt( s2-1/P^2*(s1.^2) );
y = fg./denom;
function P = iter_filtre(P,f, nb_iter)

[N,NN] = size(P);

for i=1:nb_iter
	P(1:N,i+1) = ifft( fft(f).*fft(P(1:N,i)) );
end;

P(N+1,1:nb_iter+1) = P(1,1:nb_iter+1);
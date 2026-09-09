% Multiplication de grands entiers.
function res = mult_entier(x,y,b)
N = length(x);
% ajout de zéros pour convolution acyclique
x = [x;zeros(N,1)];
y = [y;zeros(N,1)];
% calcule la convolution 
res = round( real( ifft(fft(x).*fft(y)) ) );
% enlève les retenues
for i=1:2*N-1
	q = floor(res(i)/b);
	res(i) = res(i)-q*b;
	res(i+1) = res(i+1)+q;
end
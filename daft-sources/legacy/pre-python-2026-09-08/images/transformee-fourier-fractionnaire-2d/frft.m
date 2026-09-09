function res = frft(f, alpha)

N = length(f);
w = exp( -2i*pi*alpha/N );

res = czt(f,N,w,1);
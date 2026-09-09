function fd = com_fd(x,y)

[N,NN] = size(x);

fd = complex(x,y);
fd = fft(fd);
% on suprime le 1er élément 
% (pour rendre le 1er fd invariant par translation)
fd = fd(2:N);
% on normalise 
% (pour rendre le fd invariant par similitude)
fd = fd/fd(1);
fd = fd(2:N-1);
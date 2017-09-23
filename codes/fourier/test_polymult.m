%%
% test for matrix multiplication

a = [1 0 1]';  % 1+X^2
b = [0 2 0 1]'; % 2X+X^4

% target size to avoid periodization
d = length(a)+length(b)-1;

A = [a;zeros(d-length(a),1)];
B = [b;zeros(d-length(b),1)];

% should be (1+X^2)*(2X+X^4) = 1+2X+2X^3+X^5
ifft(fft(A).*fft(B))
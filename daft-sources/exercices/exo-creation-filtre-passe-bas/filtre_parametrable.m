function f = filtre_parametrable(N,eps)
P1 = floor(eps*N/4); P = 2*P1+1; % P doit etre impair
trans = [1]; if(P~=1) trans = (cos((0:P-1)'*pi/(P-1))+1)/2; end;
f = [ones(N/4-P1,1); trans; zeros(N/2-P,1); trans(P:-1:1); ones(N/4-P1-1,1)];
f = real(ifft(f));
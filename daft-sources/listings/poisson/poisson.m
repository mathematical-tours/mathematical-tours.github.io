% quelques constantes
N = 30; h = 1/N; nb_iter = 30;
M = zeros(N+1, N+1);
% calcule du membre de droite
f_val = zeros(N-1, N-1);
% on commence avec x=h (seulement les points du centre)
for i=1:N-1	
    for j=1:N-1
        x = i*h; y = j*h;
        f_val(i,j) = f(x,y);
    end
end
% ajout des termes de bord
for i=1:N-1
    x = i*h;
    f_val(i,1)   = f_val(i,1)   - 1/h^2 * f_0y(x);
    f_val(i,N-1) = f_val(i,N-1) - 1/h^2 * f_1y(x);
    f_val(1,i)   = f_val(1,i)   - 1/h^2 * f_x0(x);
    f_val(N-1,i) = f_val(N-1,i) - 1/h^2 * f_x1(x);    
end
% on rend la matrice impaire
ff = [zeros(N-1,1),f_val,zeros(N-1,1),-f_val(:,N-1:-1:1)];
ff = [zeros(1,2*N); ff; zeros(1,2*N); -ff(N-1:-1:1,:)];
% on calcule la FFT
ff = fft2(ff);
% on calcule les vecteurs propres :
d = -4/h^2 * sin( (0:2*N-1) * pi/(2*N)).^2;
% on résout le système d*u+u*d = ff
for i=1:2*N
    for j=1:2*N
        s = d(i) + d(j);
        if s==0
            s=1;	% éviter la division par 0
        end
        ff(i,j) = ff(i,j) / s;
    end
end
% on calcule la transformée inverse
ff = real( ifft2( ff ) );
% on extrait la solution 
u = zeros(N+1,N+1);
u(2:N, 2:N) = ff(2:N,2:N);
% on remet les termes du bord
for i=1:N+1
    x = (i-1)*h;
    u(i,1) = f_0y(x);
    u(i,N+1) = f_1y(x);
    u(1,i) = f_x0(x);
    u(N+1,i) = f_x1(x);    
end  
surf(u);
title('Résolution par FFT');
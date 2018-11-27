function Fr = ComputeRadialPofile(F)

n = size(F,1);

% compute radial profile
t = -n/2:n/2-1;
[X,Y] = meshgrid(t,t);
nr = n/4; nt = n/4;
[R,T] = meshgrid(linspace(0,.98*n/2,nr), linspace(0,2*pi,nt));
Xi = R.*cos(T); Yi = R.*sin(T);
U = interp2(X,Y,fftshift(F),Xi(:),Yi(:));
U = reshape(U,size(Xi));
Fr = mean(U,1);
Fr = Fr(11:end);
Fr = Fr(:);

end
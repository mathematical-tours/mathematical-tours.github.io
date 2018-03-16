function Y = sphrharm(phi,th,L,M)

% phi in [-pi/2,pi/2]
% th in [-pi,pi]

% PHI  =linspace(0,  pi,RES(2));  
% THETA=linspace(0,2*pi,RES(1));  % Azimuthal/Longitude/Circumferential


% [Xm,Ym,Zm]=sph2cart(THETA,PHI-pi/2,abs(Ymn).^2);
% [Xr,Yr,Zr]=sph2cart(THETA,PHI-pi/2,real(Ymn).^2);
% [Xi,Yi,Zi]=sph2cart(THETA,PHI-pi/2,imag(Ymn).^2);

% th = th + pi;
phi = phi + pi/2;

Lmn=legendre(L,cos(phi));
if L~=0
    Lmn=squeeze(Lmn(abs(M)+1,:,:));
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Y = C*Lmn.*exp(1i*M*th);


if M>=0
    Y = real(Y);
else
    Y = imag(Y);
end


Y = Y/norm(Y);

end
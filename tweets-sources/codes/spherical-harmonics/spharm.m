% This function generates the Spherical Harmonics basis functions of degree
% L and order M.
%
% SYNTAX: [Ymn,THETA,PHI,X,Y,Z]=spharm4(L,M,RES,PLOT_FLAG);
%
% INPUTS:
%
% L         - Spherical harmonic degree, [1x1]
% M         - Spherical harmonic order,  [1x1]
% RES       - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1] 
% PLOT_FLAG - Binary flag to generates a figure of the spherical harmonic surfaces (DEFAULT=1)
%
%
% OUTPUTS:
% 
% Ymn   - Spherical harmonics coordinates, [RES(1) x RES(2)]
% THETA - Circumferential coordinates, [RES(1) x RES(2)]
% PHI   - Latitudinal coordinates, [RES(1) x RES(2)]
% X,Y,Z - Cartesian coordinates of magnitude, squared, spherical harmonic surface points, [RES(1) x RES(2)]
%
%
% NOTE: It is very important to keep the various usages of THETA and PHI
% straight.  For this function THETA is the Azimuthal/Longitude/Circumferential 
% coordinate and is defined on the interval [0,2*pi], whereas PHI is the 
% Altitude/Latitude/Elevation and is defined on the interval [0,pi].  Also note that 
% the conversion to cartesian coordinates requires that PHI be offset by pi/2 so 
% that the conversion is on the interval [-pi/2,pi/2].
%
% DBE 2005/09/30

function [Ymn,THETA,PHI,Xm,Ym,Zm]=spharm4(L,M,RES,PLOT_FLAG);

% Define constants (REQUIRED THAT L(DEGREE)>=M(ORDER))
if nargin==0
  L=3;   % DEGREE
  M=2;   % ORDER
  RES=[55 55];
end

if nargin<3
  RES=[25 25];  
  PLOT_FLAG=1;
end

if nargin<4
  PLOT_FLAG=1;
end

if L<M, error('The ORDER (M) must be less than or eqaul to the DEGREE(L).'); end

THETA=linspace(0,2*pi,RES(1));  % Azimuthal/Longitude/Circumferential
PHI  =linspace(0,  pi,RES(2));  % Altitude /Latitude /Elevation
[THETA,PHI]=meshgrid(THETA,PHI);

Lmn=legendre(L,cos(PHI));

if L~=0
  Lmn=squeeze(Lmn(M+1,:,:));
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Ymn=C*Lmn.*exp(i*M*THETA);

[Xm,Ym,Zm]=sph2cart(THETA,PHI-pi/2,abs(Ymn).^2);
[Xr,Yr,Zr]=sph2cart(THETA,PHI-pi/2,real(Ymn).^2);
[Xi,Yi,Zi]=sph2cart(THETA,PHI-pi/2,imag(Ymn).^2);
% [Xp,Yp,Zp]=sph2cart(THETA,PHI-pi/2,angle(Ymn).^2);

if PLOT_FLAG
f=figure; axis off; hold on;
  axes('position',[0.0500 0 0.2666 1]); 
    surf(Xm,Ym,Zm); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0.3666 0 0.2666 1]); 
    surf(Xr,Yr,Zr); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0.6833 0 0.2666 1]); 
    surf(Xi,Yi,Zi); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0 0.9 1 0.1]); axis off;
    t(1)=text(0.50,0.25,'Spherical Harmonics','HorizontalAlignment','Center');
  axes('position',[0 0 1 0.1]); axis off;
    t(2)=text(0.20,0.25,['|Y^',num2str(M),'_',num2str(L),'|^2'],'HorizontalAlignment','Center');
    t(3)=text(0.50,0.25,['Real(Y^',num2str(M),'_',num2str(L),')^2'],'HorizontalAlignment','Center');
    t(4)=text(0.80,0.25,['Imag(Y^',num2str(M),'_',num2str(L),')^2'],'HorizontalAlignment','Center');
  setfig(gcf,10,5,12);
end

return
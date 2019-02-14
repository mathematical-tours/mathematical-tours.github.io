function [g,N] = GPA(f, sigmar, W, eps, flag)
% Gaussian Bilateral filter:
% [g,Ng] = GPA(f, sigmar, sigmas, eps, 'Gauss')
% f             : input image 
% sigmar        : width of range Gaussian
% sigmas        : width of spatial Gaussian
% eps           : kernel approximation accuracy
% g             : output image
% Ng            : approximation order
%
% Box bilateral filter:
% [b,Nb] = GPA(f, sigmar, B, eps, 'box')
% f             : input image 
% sigmar        : width of range Gaussian
% B             : width of box kernel
% eps           : kernel approximation accuracy
% g             : output image
% Nb            : approximation order
%
if strcmp(flag,'Gauss')
    L=round(3*W);
	Hs=fspecial('gaussian',2*L+1,W);
elseif  strcmp(flag,'box')
    L=W;
	Hs=fspecial('average',2*L+1);
else
	error('not enough arguments');
end 
T = 128;
N = estN(sigmar,T,eps);
f=padarray(f,[L,L]);
H=(f-T)/sigmar;     
F=exp(-0.5*H.^2);   
G=ones(size(H));
P=zeros(size(H));  
Q=zeros(size(H));   
Fbar=imfilter(F,Hs);   
for n = 1 : N
	Q=Q+G.*Fbar;
    F=H.*F/sqrt(n);
    Fbar=imfilter(F,Hs);
    P=P+G.*Fbar*sqrt(n);
    G=H.*G/sqrt(n);
end
g= T+sigmar*(P(L+1:end-L,L+1:end-L) ...
    ./Q(L+1:end-L,L+1:end-L));      
g(g<0)=0;
g(g>255)=255;
end
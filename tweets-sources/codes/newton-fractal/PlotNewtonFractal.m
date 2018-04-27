function [C1,x,y] = PlotNewtonFractal(Phi,R,N,A,u)

% Phi : map to iterate
% R : roots
% A: width of window
% u: center

% colormap
cm = [[1 0 0]; [0 1 0]; [0 0 1]; ...
    [0 1 1]; [1 0 1]; [1 1 0]
    [1/2 1 1]; [1 1/2 1]; [1 1 1/2]
    [1/2 1 1/2]; [1/2 1/2 1]; [1 1/2 1/2]]';

x = linspace(real(u)-A,real(u)+A,N);
y = linspace(imag(u)-A,imag(u)+A,N);
% grid generation
[Y,X] = meshgrid(y,x,N);
Z = X + 1i*Y;

niter = 150;
for i=1:niter
    progressbar(i,niter);
    Z = Phi(Z);
end
% distance to roots
D = zeros(N,N,length(R));
for k=1:length(R)
    D(:,:,k) = abs(Z-R(k));
end
[D1,C] = min(D,[],3);
% colormap
C1 = zeros(N,N,3);
for k=1:length(R)
    for i=1:3
        C1(:,:,i) = C1(:,:,i) + (C==k) * cm(i,k);
    end
end

end
%%
% Displays the bassin of attraction of Newton method on the complex plane.

rep = 'results/';
[~,~] = mkdir(rep);

if not(exist('q'))
    q = 3;
end
if not(exist('theta'))
    theta = 1.4;
    rho = 1;
end

f  = @(z)z.^q-1;
Df = @(z)q*z.^(q-1);
a = 1.5;
R = exp(2i*pi/q*(0:q-1)');

tau = 0;
Newton = @(z)z - ( 1 - rho*exp(1i*theta) ) * f(z)./Df(z);
Phi = @(z)(1-tau)*Newton(z)+tau*z;

% grid
N = 512;
x = linspace(-a,a,N);
[Y,X] = meshgrid(x,x,N);
Z = X + 1i*Y;

niter = 400;
for i=1:niter
    progressbar(i,niter);
    Z = Phi(Z);
end

% roots
D = zeros(N,N,length(R));
for k=1:length(R)
    D(:,:,k) = abs(Z-R(k));
end
[D1,C] = min(D,[],3);

clf;
imagesc(C); 
axis image; axis off;

% colormap 
cm = [[1 0 0]; [0 1 0]; [0 0 1]; ...
      [0 1 1]; [1 0 1]; [1 1 0]
      [1/2 1 1]; [1 1/2 1]; [1 1 1/2]
      [1/2 1 1/2]; [1/2 1/2 1]; [1 1/2 1/2]]';
C1 = zeros(N,N,3);
for k=1:length(R)
    for i=1:3
        C1(:,:,i) = C1(:,:,i) + (C==k) * cm(i,k);
    end
end
% imwrite(C1, [rep 'newton-fractal-' num2str(q) '.png'], 'png');
imwrite(C1, [rep 'newton-fractal-' num2str(q) '-rho' num2str(rho) '.png'], 'png');



clf; hold on;
imagesc(x,x,permute(C1, [2 1 3])); 
axis image; axis off;
plot(R,'k.', 'MarkerSize', 40);
drawnow;
% saveas(gcf,  [rep 'newton-fractal-' num2str(q) '-pt.png'], 'png');
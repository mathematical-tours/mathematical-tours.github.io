%%
% Brenier in heat

n = 1024/2;
H = @(u)u;
H = @(u)atan(u/200);

psi = (randn(n,1)+1/2)*2;
t = (0:n-1)'/n;
u0 = (abs(t-.5)<.25); 
u0 = cos(t*2*pi);

b = 10;
a = -b;
u0 = (u0-min(u0)) / (max(u0)-min(u0));
u0 = (b-a)*u0+a;


psi0 = 0.5*cos(2*pi*(t-0.25))+0.2*cos(10.*pi*(t-0.2));
% theta0 =1.5*cos(2*pi*(t-0.25))+0.2*cos(10.*pi*(t-0.2));
psi0 = psi0 * 8;
Delta = @(f)(2*f - f([2:end 1]) - f([end 1:end-1]) )/2;

T = .1;
h = .05/n^2;
h = .003/n^2;
niter = round(T/h);
niter = 5000;
ndisp = max(1,ceil(niter/20));

clf; hold on;
psi = psi0;
Psi = [];
for i=1:niter
    Psi(:,end+1) = psi;
    if mod(i,ndisp)==1
        c = (i-1)/(niter-1);
        plot(psi, 'LineWidth', 2, 'Color', [1-c,0,c]);
        drawnow;
    end
    psi = psi - Delta(psi).*cos( pi*psi ).^2;
end
axis tight;

%% 
% resample for display
D = 
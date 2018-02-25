%%
% Brenier in heat

n = 1024/4;
H = @(u)u;
H = @(u)atan(u/200);

u = (randn(n,1)+1/2)*2;
t = (0:n-1)'/n;
u0 = (abs(t-.5)<.25); 
u0 = cos(t*2*pi);

b = 10;
a = -b;
u0 = (u0-min(u0)) / (max(u0)-min(u0));
u0 = (b-a)*u0+a;


theta0 = 0.5*cos(2*pi*(t-0.25))+0.2*cos(10.*pi*(t-0.2));
theta0 =1.5*cos(2*pi*(t-0.25))+0.2*cos(10.*pi*(t-0.2));
%  pi/2*theta=Arctg u
u0 = tan(2/pi*theta0);

Delta = @(f)(2*f - f([2:end 1]) - f([end 1:end-1]) ) * 2*n^2;

T = .1;
h = .05/n^2;
h = .003/n^2;
niter = round(T/h);
ndisp = max(1,ceil(niter/10));

clf; hold on;
u = u0;
for i=1:niter
    if mod(i,ndisp)==1
        c = (i-1)/(niter-1);
        theta = atan(u)*2/pi;
        plot(theta/max(abs(theta)), 'LineWidth', 2, 'Color', [1-c,0,c]);
        %plot(u, 'LineWidth', 2, 'Color', [1-c,0,c]);
        drawnow;
    end
    u = u - Delta(H(u))*h;
end
axis tight;
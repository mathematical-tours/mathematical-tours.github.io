%%
% Contrast iterative mean vs iterative median.

% # nn
k = 3;
% size of signal
n = 512*2;

x = 1:n;
s = -k:k;
[S,X] = meshgrid(s,x);
I = mod(X + S - 1,n)+1;

t = (0:n-1)'/n;
sigma = 1;
f0 = sin(t*2*pi) + sigma*randn(n,1);

oper = @(u)median(u,2);
oper = @(u)mean(u,2);

niter = 500;
nbdisp = 50;
ndisp = round( linspace(1,sqrt(niter),nbdisp).^2 );
ndisp = 1:niter;

q = max(abs(f0));
f = f0;
a = 1;
for it=1:max(ndisp);
    f = oper(f(I));
    if it<=length(ndisp) && it==ndisp(a)
        c = (a-1)/(nbdisp-1);
        clf;
        plot(f, 'LineWidth', 2, 'Color', [c 0 1-c]);
        axis([1 n -q q]); axis off;
        drawnow;
        a = a+1;
    end
end
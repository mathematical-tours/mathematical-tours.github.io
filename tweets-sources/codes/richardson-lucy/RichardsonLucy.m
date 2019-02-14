%%
% Test for mirror-descent deconvolution under positivity constraints.

n = 256*10; 

% load filter
t = [0:n/2,-n/2+1:-1]'/n;

% Gaussian kernel
s = .05; 
k = exp(-t.^2/(2*s^2)); k = k/sum(k);
k = k/sum(k);
% ED kernel
s = .05; 
k = -abs(t); 

% be sure it is positive
kf = max(fft(k),0);
% convolution operator
K = @(x)real(ifft(kf.*fft(x)));

% data
x0 = double(abs(t-1/2)<.08);
x0 = x0/sum(x0);

% energy <h*(x-x0),x-x0>
dotp = @(x,y)sum(x(:).*y(:));
E = @(x)dotp(K(x-x0),x-x0)/2;

niter = 100000;
q = 80; % #display
ndisp = round(linspace(1,niter,q));
kdisp = 1;


tau = 50; %% NEEDS TO BE TUNED
x = ones(n,1)/n; 
for it=1:niter
    G = K(x-x0);
    x = x .* exp( -tau*G );
    x = x/sum(x);
    Elist(it)=E(x);
    if it==ndisp(kdisp)
        clf; hold on;
        plot(x0, 'r', 'LineWidth', 2); plot(x0, 'r');
        plot(x, 'k', 'LineWidth', 2); axis([1 n 0 max(x0)*1.05]); box on;
        drawnow;
        kdisp = kdisp+1;
    end
end
return;
clf;
plot(log(Elist), 'LineWidth', 2); 
axis tight;

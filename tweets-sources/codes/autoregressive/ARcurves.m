%% 
% Tests for 2D AR processes
% X_{t+1} = a*X_t + b*X_{t-1} + W_t

% impose : 
% a = 2*r;
% a = -a^2 / 4;


% width^2 of the Gaussian-like covariance 
sigma = 100^2; 

[a,Sigma,rlist] = fit_ar2(sigma);
b = a(2); a = a(1); 
gamma = @(k,r)r.^abs(k) .* ( 1+(1-r^2)/(1+r^2)*abs(k) );
t = 1:1000;
clf; plot(t,gamma(t,a/2)); % display correlation


crand = @(a,b)randn(a,b)+1i*rand(a,b);
n = 10000;
p = 10; % #particles
wu = 5000; % warmup stage 

x = zeros(2,p);
for i=1:n+wu
    x(end+1,:) = a*x(end,:) + b*x(end-1,:) + crand(1,p);
end
x = x(wu:end,:);
plot(x, '-')
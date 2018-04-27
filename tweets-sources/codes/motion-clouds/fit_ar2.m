function [a,Sigma,rlist] = fit_ar2(sigma)

% fit_ar2 - fit the parameter of an AR(2) process to match some given covariance width
%
%   a = fit_ar2(sigma);
%
%   Fit the parameter (a,b) of the AR(2) process
%       X_{t+1} = a(1)*X_t + a(2)*X_{t-1} + W_t
%   so that the width of the covariance is near to sigma.
%
%   We impose the relation
%       a(2) = -a(1)^2 / 4   
%
%   Copyright (c) 2013 Gabriel Peyre

% covariance as function of r=a(2)/2.
gamma = @(k,r)r.^abs(k) .* ( 1+(1-r^2)/(1+r^2)*abs(k) );

klist = -5000:5000;
Q = 400;
rlist = linspace(.01, .99, Q);

Sigma = zeros(Q,1);
for i=1:Q
    r = rlist(i);
    g = gamma(klist,r);
    Sigma(i) = sqrt( sum( g .* klist.^2 ) );
end

[~,i] = min(abs(Sigma-sigma));
r = rlist(i);

a(1) = 2*r;
a(2) = -a(1)^2 / 4;

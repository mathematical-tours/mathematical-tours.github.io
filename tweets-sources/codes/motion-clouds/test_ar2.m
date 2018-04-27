r = .5;
gamma = @(k,r)r.^abs(k) .* ( 1+(1-r^2)/(1+r^2)*abs(k) );

k = -20:20;
plot(k, gamma(k,r) );

k = -5000:5000;
r = linspace(.02, .98,200);

sigma = [];
for i=1:length(r)
    g = gamma(k,r(i));
    sigma(i) = sqrt( sum( g .* k.^2 ) );
end

plot(r, sigma);

plot(r, log(sigma));
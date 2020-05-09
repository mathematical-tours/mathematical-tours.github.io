%% 
% differentiate along a curve

curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );

N = 1024;
c0 = exp(2i*pi*(0:N-1)'/N);
c0 = c0 ./ max(abs(real(c0)), abs(imag(c0)));
n = 50; % #zeros
c = resample(c0, n);

P = poly(c);
nder = 25;
P1 = P;
for k=1:nder
    P1 = polyder(P1);
end
cD = roots(P1);

clf; hold on;
plot(c([1:end 1]), 'b'); 
plot(cD, 'r.', 'MarkerSize', 10); 
axis equal; axis([-1 1 -1 1]*1.05);
axis off;


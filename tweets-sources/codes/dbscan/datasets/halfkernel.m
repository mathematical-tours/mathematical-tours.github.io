function data = halfkernel(N, minx, r1, r2, noise, ratio)

if nargin < 1
    N = 1000;
end
if mod(N,2) ~= 0
    N = N + 1;
end
if nargin < 2
    minx = -20;
end
if nargin < 3
    r1 = 20;
end
if nargin < 4
    r2 = 35;
end
if nargin < 5
    noise = 4;
end
if nargin < 6
    ratio = 0.6;
end

phi1 = rand(N/2,1) * pi;
inner = [minx + r1 * sin(phi1) - .5 * noise  + noise * rand(N/2,1) r1 * ratio * cos(phi1) - .5 * noise + noise * rand(N/2,1) ones(N/2,1)];
    
phi2 = rand(N/2,1) * pi;
outer = [minx + r2 * sin(phi2) - .5 * noise  + noise * rand(N/2,1) r2 * ratio * cos(phi2) - .5 * noise  + noise * rand(N/2,1) zeros(N/2,1)];

data = [inner; outer];
end
function data = twospirals(N, degrees, start, noise,rep)
% Generate "two spirals" dataset with N*rep instances.
% degrees controls the length of the spirals
% start determines how far from the origin the spirals start, in degrees
% noise displaces the instances from the spiral. 
%  0 is no noise, at 1 the spirals will start overlapping

    if nargin < 1
        N = 200;
    end
    if nargin < 2
        degrees = 570;
    end
    if nargin < 3
        start = 90;
    end
    if nargin < 5
        noise = 0.2;
    end  
    
    deg2rad = (2*pi)/360;
    start = start * deg2rad;

    N1 = floor(N/2);
    N2 = N-N1;
    
    n = start + sqrt(rand(N1,1)) * degrees * deg2rad;
    d1=[];
    for i=1:rep
        d1 = [d1;[-cos(n).*n + rand(N1,1)*noise sin(n).*n+rand(N1,1)*noise zeros(N1,1)]];
    end
    
    n = start + sqrt(rand(N1,1)) * degrees * deg2rad;    
    d2=[];
    for i=1:rep
        d2 = [d2;[cos(n).*n+rand(N2,1)*noise -sin(n).*n+rand(N2,1)*noise ones(N2,1)]];
    end
    
    data = [d1;d2];
end
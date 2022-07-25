function data = outlier(N, r, dist, outliers, noise)

    if nargin < 1
        N = 600;
    end
    if nargin < 2
        r = 20;
    end
    if nargin < 3
        dist = 30;
    end
    if nargin < 4
        outliers = 0.04;
    end
    if nargin < 5
        noise = 5;
    end

    N1 = round(N * (.5-outliers));
    N2 = N1;
    N3 = round(N * outliers);
    N4 = N-N1-N2-N3;

    phi1 = rand(N1,1) * pi;
    r1 = sqrt(rand(N1,1))*r;
    P1 = [-dist + r1.*sin(phi1) r1.*cos(phi1) zeros(N1,1)];

    phi2 = rand(N2,1) * pi;
    r2 = sqrt(rand(N2,1))*r;
    P2 = [dist - r2.*sin(phi2) r2.*cos(phi2) 3*ones(N2,1)];    
    
    P3 = [rand(N3,1)*noise dist+rand(N3,1)*noise 2*ones(N3,1)];    
    
    P4 = [rand(N4,1)*noise -dist+rand(N4,1)*noise ones(N4,1)];
    
    data = [P1; P2; P3; P4];

end
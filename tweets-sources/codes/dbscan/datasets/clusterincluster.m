function data = clusterincluster(N, r1, r2, w1, w2, arms)

    if nargin < 1
        N = 1000;
    end
    if nargin < 2
        r1 = 1;
    end
    if nargin < 3
        r2 = 5*r1;
    end
    if nargin < 4
        w1 = 0.8;
    end
    if nargin < 5
        w2 = 1/3;
    end
    if nargin < 6
        arms = 64;
    end
    
    data = [];
    
    N1 = floor(N/2);
    N2 = N-N1;
    
    phi1 = rand(N1,1) * 2 * pi;
    dist1 = r1 + randint(N1,1,3)/3 * r1 * w1;
    d1 = [dist1 .* cos(phi1) dist1 .* sin(phi1) zeros(N1,1)];

    perarm = round(N2/arms);
    N2 = perarm * arms;
    radperarm = (2*pi)/arms;
    phi2 = ((1:N2) - mod(1:N2, perarm))/perarm * (radperarm);
    phi2 = phi2';
    dist2 = r2 * (1 - w2/2) + r2 * w2 * mod(1:N2, perarm)'/perarm;
    d2 = [dist2 .* cos(phi2) dist2 .* sin(phi2) ones(N2,1)];    
    
    data = [d1;d2];   

    scatter(data(:,1), data(:,2), 20, data(:,3)); axis square;
end
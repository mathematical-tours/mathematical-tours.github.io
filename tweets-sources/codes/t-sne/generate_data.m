function [X, labels, t] = generate_data(dataname, n, noise)
%GENERATE_DATA Generates an artificial dataset (manifold)
%
%	[X, labels, t] = generate_data(dataname, n, noise)
%
% Generates an artificial dataset. Possible datasets are: 'swiss' for the Swiss roll
% dataset, 'helix' for the helix dataset, 'twinpeaks' for the twinpeaks dataset,
% '3d_clusters' for the 3D clusters dataset, and 'intersect' for the intersecting
% dataset. The variable n indicates the number of datapoints to generate
% (default = 1000). The variable noise indicates the amount of noise that
% is added to the data (default = 0.05). The function returns the
% high-dimensional dataset in X, and corresponding labels in labels. In
% addition, the function returns the coordinates of the datapoints on the
% underlying manifold in t.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology



if ~exist('n', 'var')
    n = 1000;
end
if ~exist('noise', 'var')
    noise = 0.05;
end

switch dataname
    case 'swiss'
        t = (3 * pi / 2) * (1 + 2 * rand(n, 1));
        height = 30 * rand(n, 1);
        X = [t .* cos(t) height t .* sin(t)] + noise * randn(n, 3);
        %labels = uint8(t);
        labels = rem(sum([round(t / 2) round(height / 12)], 2), 2);
        t = [t height];
        
    case 'brokenswiss'
        t = [(3 * pi / 2) * (1 + 2 * rand(ceil(n / 2), 1) * .4); (3 * pi / 2) * (1 + 2 * (rand(floor(n / 2), 1) * .4 + .6))];
        height = 30 * rand(n, 1);
        X = [t .* cos(t) height t .* sin(t)] + noise * randn(n, 3);
        labels = uint8(t);
        %labels = rem(sum([round(t / 2) round(height / 12)], 2), 2);
        t = [t height];
        
    case 'changing_swiss'
        r = zeros(1, n);
        for i=1:n
            pass = 0;
            while ~pass
                rr = rand(1);
                if rand(1) > rr
                    r(i) = rr;
                    pass = 1;
                end
            end
        end
        t = (3 * pi / 2) * (1 + 2 * r);
        height = 21 * rand(1, n);
        X = [t .* cos(t); height; t .* sin(t)]' + noise * randn(n, 3);
        %labels = uint8(t)';
        labels = rem(sum([round(t / 2); round(height / 10)], 1), 2)';
        
    case 'helix'
        t = [1:n]' / n;
        t = t .^ (1.0) * 2 * pi;
        X = [(2 + cos(8 * t)) .* cos(t) (2 + cos(8 * t)) .* sin(t) sin(8 * t)] + noise * randn(n, 3);
        %labels = uint8(t);
        labels = rem(round(t * 1.5), 2);
        
    case 'twinpeaks'
        inc = 1.5 / sqrt(n);
        [xx2, yy2] = meshgrid(-1:inc:1);
        xy = 1 - 2 * rand(2, n);
        X = [xy; sin(pi * xy(1,:)) .* tanh(3 * xy(2,:))]' + noise * randn(n, 3);
        X(:,3) = X(:,3) * 10;
        t = xy';
        %labels = uint8(X(:,3));
        labels = rem(sum(round((X + repmat(min(X, [], 1), [size(X, 1) 1])) ./ 10), 2), 2);
        
    case '3d_clusters'
        numClusters = 5;
        centers = 10 * rand(numClusters, 3);
        D = L2_distance(centers', centers');
        minDistance = min(D(D > 0));
        k = 1;
        n2 = n - (numClusters - 1) * 9;
        X = repmat(0, [n 3]);
        labels = repmat(0, [n 1]);
        for i=1:numClusters
            for j=1:ceil(n2 / numClusters)
                X(k, 1:3) = centers(i, 1:3) + (rand(1, 3) - 0.5) * minDistance / sqrt(12);
                labels(k) = i;
                k = k + 1;
            end
        end
        X = X + noise * randn(size(X, 1), 3);
        t = [];
        
    case 'intersect'
        t = [1:n]' ./ n .* (2 * pi);
        x = cos(t);
        y = sin(t);
        height = rand(length(x), 1) * 5;
        X = [x x .* y height] + noise * randn(n, 3);
        %labels = uint8(5 * t);
        labels = rem(sum([round(t / 2) round(height / 2)], 2), 2);
        
    case 'difficult'
        % Generate underlying manifold
        no_dims = 5;
        no_points_per_dim = round(n ^ (1 / no_dims));
        l = linspace(0, 1, no_points_per_dim);
        t = combn(l, no_dims);
        
        % Generate high-dimensional dataset
        X = [cos(t(:,1)) tanh(3 * t(:,2)) t(:,1) + t(:,3) t(:,4) .* sin(t(:,2)) sin(t(:,1) + t(:,5)) t(:,5) .* cos(t(:,2)) t(:,5) + t(:,4) t(:,2) t(:,3) .* t(:,4) t(:,1)];
        X = X + noise * randn(size(X));
        
        % Generate labels for dataset (2x2x2x2x2 checkerboard pattern)
        tt = 1 + round(t);
        labels = rem(sum(tt, 2), 2);
        
    otherwise
        error('Unknown dataset name.');
end
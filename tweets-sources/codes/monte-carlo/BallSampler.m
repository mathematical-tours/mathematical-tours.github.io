function X = BallSampler(d,n)

% draw on sphere
X = randn(d,n); 
X = X./repmat( sqrt(sum(X.^2,1)), [d 1] );
% map to ball
X = X .* repmat( rand(1,n).^(1/d) , [d 1] );

end
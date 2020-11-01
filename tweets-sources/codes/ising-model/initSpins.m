function spin = initSpins(numSpinsPerDim, p)
%INITSPINS Initialize a configuration of spins.
%   spin = INITSPINS(numSpinsPerDim, p) returns a configuration of spins
%   with |numSpinsPerDim| spins along each dimension and a proportion |p|
%   of them pointing upwards. |spin| is a matrix of +/- 1's.
%   Copyright 2017 The MathWorks, Inc.
spin = sign(p - rand(numSpinsPerDim, numSpinsPerDim));

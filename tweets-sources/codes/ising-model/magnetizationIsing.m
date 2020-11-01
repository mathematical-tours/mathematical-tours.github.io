function Mmean = magnetizationIsing(spin)
%MAGNETIZATIONISING Magnetization of a configuration of spins.
%   Mmean = MAGNETIZATIONISING(spin) returns the magnetization of the
%   configuration |spin|. |spin| is a matrix of +/- 1's.
%   Copyright 2017 The MathWorks, Inc.
Mmean = mean(spin(:));
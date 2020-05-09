function out = decomposeATriangle(in)
%decomposeATriangle Decompose type A triangle.
%
%   t = decomposeATriangle(in) decomposes a type A triangle into a type A
%   and a type B-prime triangle. The input is a one-row table as returned
%   by aTriangle, and the output is a two-row table.

%   Copyright 2018 The MathWorks, Inc.

out = [ ...
    aTriangle(in.Left,in.Right,[])
    bpTriangle([],in.Apex,in.Left) ];

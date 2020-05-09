function out = decomposeApTriangle(in)
%decomposeApTriangle Decompose type A-prime triangle.
%
%   t = decomposeApTriangle(in) decomposes a type A-prime triangle into a
%   type A-prime and a type B triangle. The input is a one-row table as
%   returned by apTriangle, and the output is a two-row table.

%   Copyright 2018 The MathWorks, Inc.

out = [ ...
    apTriangle(in.Right,[],in.Left)
    bTriangle([],in.Right,in.Apex) ];

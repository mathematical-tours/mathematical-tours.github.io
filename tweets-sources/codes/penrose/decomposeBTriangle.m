function out = decomposeBTriangle(in)
%decomposeBTriangle Decompose type B triangle.
%
%   t = decomposeBTriangle(in) decomposes a type B triangle into a type A
%   triangle, a type B triangle, and a type B-prime triangle. The input is
%   a one-row table as returned by bTriangle, and the output is a three-row
%   table.

%   Copyright 2018 The MathWorks, Inc.

t_b = bTriangle([],in.Right,in.Apex);
t_a = aTriangle(t_b.Apex,in.Apex,[]);
t_bp = bpTriangle(t_a.Right,in.Left,t_a.Apex);

out = [t_b ; t_a ; t_bp];


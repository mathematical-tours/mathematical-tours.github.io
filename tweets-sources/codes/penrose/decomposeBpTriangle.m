function out = decomposeBpTriangle(in)
%decomposeBpTriangle Decompose type B-prime triangle.
%
%   t = decomposeBpTriangle(in) decomposes a type B-prime triangle into a
%   type B-prime triangle, a type A-prime triangle, and a type B triangle.
%   The input is a one-row table as returned by bpTriangle, and the output
%   is a three-row table.

%   Copyright 2018 The MathWorks, Inc.

t_bp = bpTriangle([],in.Apex,in.Left);
t_ap = apTriangle(t_bp.Apex,[],in.Apex);
t_b = bTriangle(t_ap.Left,t_ap.Apex,in.Right);

out = [t_bp ; t_ap ; t_b];



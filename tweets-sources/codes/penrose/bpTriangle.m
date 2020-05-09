function t = bpTriangle(Apex,Left,Right)
%bpTriangle Representation of type B-prime triangle in table form.
%   t = bpTriangle(apex,left,right) returns a type B-prime triangle
%   represented as a one-row table. The triangle vertices (apex, left, and
%   right) are complex numbers. If one of the triangle vertices is empty,
%   it is computed automatically.
%
%   EXAMPLE
%
%   Compute a type B-prime triangle with apex at (0,1) and left base vertex
%   at (0,0) in the complex plane.
%
%       t = bpTriangle(1i,0,[])

%   Copyright 2018 The MathWorks, Inc.

[Apex,Left,Right] = isoscelesTriangle(Apex,Left,Right,108);
Type = categorical("Bp");
t = table(Apex,Left,Right,Type);
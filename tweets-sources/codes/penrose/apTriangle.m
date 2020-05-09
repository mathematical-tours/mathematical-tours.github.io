function t = apTriangle(Apex,Left,Right)
%aTriangle Representation of type A-prime triangle in table form.
%   t = aTriangle(apex,left,right) returns a type A-prime triangle
%   represented as a one-row table. The triangle vertices (apex, left, and
%   right) are complex numbers. If one of the triangle vertices is empty,
%   it is computed automatically.
%
%   EXAMPLE
%
%   Compute a type A-prime triangle with apex at (0,1) and left base vertex
%   at (0,0) in the complex plane.
%
%       t = apTriangle(1i,0,[])

%   Copyright 2018 The MathWorks, Inc.

[Apex,Left,Right] = isoscelesTriangle(Apex,Left,Right,36);
Type = categorical("Ap");
t = table(Apex,Left,Right,Type);
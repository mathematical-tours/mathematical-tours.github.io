function t = bTriangle(Apex,Left,Right)
%bTriangle Representation of type B triangle in table form.
%   t = bTriangle(apex,left,right) returns a type B triangle represented as
%   a one-row table. The triangle vertices (apex, left, and right) are
%   complex numbers. If one of the triangle vertices is empty, it is
%   computed automatically.
%
%   EXAMPLE
%   Compute a type B triangle with apex at (0,1) and left base vertex at
%   (0,0) in the complex plane.
%
%       t = bTriangle(1i,0,[])

%   Copyright 2018 The MathWorks, Inc.

[Apex,Left,Right] = isoscelesTriangle(Apex,Left,Right,108);
Type = categorical("B");
t = table(Apex,Left,Right,Type);
function t = aTriangle(Apex,Left,Right)
%aTriangle Representation of type A triangle in table form.
%   t = aTriangle(apex,left,right) returns a type A triangle represented as
%   a one-row table. The triangle vertices (apex, left, and right) are
%   complex numbers. If one of the triangle vertices is empty, it is
%   computed automatically.
%
%   EXAMPLE
%   Compute a type A triangle with apex at (0,1) and left base vertex at
%   (0,0) in the complex plane.
%
%       t = aTriangle(1i,0,[])

%   Copyright 2018 The MathWorks, Inc.

[Apex,Left,Right] = isoscelesTriangle(Apex,Left,Right,36);
Type = categorical("A");
t = table(Apex,Left,Right,Type);
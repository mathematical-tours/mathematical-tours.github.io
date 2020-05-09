function [apex,left,right] = isoscelesTriangle(apex,left,right,theta)
%isoscelesTriangle Isosceles triangle.
%   [apex,left,right] = isoscelesTriangle(apex,left,right,theta) returns
%   the three vertices of an isosceles triangle given any two vertices and
%   the apex angle (in degrees). Triangle vertices are represented as
%   points in the complex plane. Specify the unknown input vertex as [].
%
%   If none of the input vertices is empty, then they are returned
%   unmodified.
%
%   EXAMPLE
%   Compute the isosceles triangle with apex at (0,1), left base vertex at
%   (0,0), and an apex angle of 36 degrees.
%
%       [apex,left,right] = isoscelesTriangle(1i,0,[],36)
%       v = [apex left right apex];
%       plot(real(v),imag(v),'LineWidth',2)
%       axis equal

%   Copyright 2018 The MathWorks, Inc.

if isempty(apex)
    base_to_side_ratio = sqrt(2 - 2*cos(deg2rad(theta)));
    base = right - left;
    apex = left + (abs(base)/base_to_side_ratio) * exp(1i * (angle(base) + deg2rad((180-theta)/2)));
    apex = scrub(apex);
    
elseif isempty(left)
    right_side = right - apex;
    left = apex + abs(right_side) * exp(1i * (angle(right_side) - deg2rad(theta)));
    left = scrub(left);
    
elseif isempty(right)
    left_side = left - apex;
    right = apex + abs(left_side) * exp(1i * (angle(left_side) + deg2rad(theta)));
    right = scrub(right);
end

function z = scrub(z)
% z = scrub(z) removes unsightly real or imaginary part.
if abs(real(z)) < 10*eps(abs(imag(z)))
    z = imag(z)*1i;
end
if abs(imag(z)) < 10*eps(abs(real(z)))
    z = real(z);
end

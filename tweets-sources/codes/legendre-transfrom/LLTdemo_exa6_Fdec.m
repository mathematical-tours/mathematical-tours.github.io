function Fdec=LLTdemo_exa6_Fdec(x)
% LLTdemo_exa6_Fdec	When x belongs to [c-a,c+a], computes the lower
%			part of a circle of radius a centered in c,
%			otherwise returns infinity.
%
% Usage:
%  Fdec=LLTdemo_exa6_Fdec(x)
%
% Input parameters:
%  x - real number
%
%  Global parameter:
%  a - radius of the circle
%  c - center of the circle
%
% Output parameter:
%  Fdec - Lower part of the circle or infinity. This is a real
%	number.
%

global a c;
if abs(x-c)<=a,
	Fdec=-sqrt(a^2-(x-c)^2);
else
	Fdec=inf;
end;


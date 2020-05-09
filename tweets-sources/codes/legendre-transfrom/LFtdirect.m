function [S,conjuguee] = LFtdirect(F,a,b,n,c,d,m)
% LFtdirect	Compute numerically the convex conjugate of a
%		univariate function given as a separate file.
%		It evaluate the function F.m at n equidistant points
%		on [a,b] and return its conjugate at m (or less)
%		equidistant points on [c,d]. It runs in quadratic
%		time.
%
% Usage:
%  [S,conjuguee] = LFtdirect(F,a,b,n,c,d,m)
%
% Input parameters:
%  F - a string containing the name of the function.
%  a,b - real number with [a,b] being the interval on which the
%	function is sampled at n points (a<b).
%  n - number of points to sampled the function.
%  c,d - real numbers with [c,d] being the interval on which we want
%	to know the conjugate.
%  m - number of points at which the conjugate is computed.
%
%
% Output parameter:
%  S - Slopes needed to define the conjugate. Removing affine parts
%	may result in needing less than m slopes. This is a column
%	vector.
%  conjuguee - numerical values of the conjugate at slopes S. This is
%	a column vector.
%
% Example:
%  LFt('f',-2,2,30,-3,3,20) % where file f.m contains the following
%			    % data:	function ff=f(x)
%			    %		 if abs(x)>2,
%			    %		   ff=(x^2-4)^2;
%			    %		 else ff=0;
%			    %		 end;
%			    %		end
%
%

X=[a:(b-a)/(n-1):b]'; % Generate the primal points
for i=1:size(X,1), Y(i)=feval(F,X(i));end; % sample the function
Y=Y';
S=[c:(d-c)/(m-1):d]'; % Generate the slopes

conjuguee = max(X * S' - Y* ones(size(S))' )'; % Compute the discrete Legendre Transform

